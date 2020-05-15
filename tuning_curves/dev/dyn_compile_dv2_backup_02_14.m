% [x, frbins, Pjoints, fr_given_as, fr_var_given_as, a_given_frs] = compile_dv2(cellid, t0s, {'lag', 0.2}, ...
%         {'frbins', 0:0.1:10}, {'dt', 0.01}, {'trialnums', []})
%
% For each timepoint in t0s, looks over trials to compile the joint
% distribution P(r,a | t=t0) where r=cellid's firing rate and a=value of accumulator
% decision variable. The firing rate on each trial is examined at time
% t0+lag. Data is obtained from cell_packager, normalized firing rates
% from Tim's DV files, and model parameters from Tim's Ephys_best_fits2.mat
%
% PARAMETERS:
% -----------
%
% cellid    An integer identifying the neuron
%
% t0s       A vector of times to examine. These will be times into the
%           model's evolution relative to alignment.
% 
% OPTIONAL PARAMETERS:
% --------------------
%
%  lag      Default 0.2; firing rates are examined this much time after the
%           correspinding time in the model and the clicks.
%
% frbins    Bins into which normalized firing rate information will be
%           divided. Default is 0:0.1:10.
%
% dt        Timestep at which to run the model. Default 0.01 (so it looks
%           nice for displaying).
%
% trialnums  Vector, a list of trials to compile over. Default is empty,
%           which means compile over all available trials.
%
% use_nans  Default is zero. If passed as 1, then firing rate rows of
%           Pjoints that did not have any info put in them will be full of
%           NaNs; if passed as 1, they will be full of zeros.
%
% alignment string matching one of align_strs in align_LUT. This will set
%           alignment of time zero of analysis (e.g., stimstart, cpokeend,
%           cpokeout)
%
% 
% RETURNS:
% --------
%
% x         Vector indexing values of accumulator decision variable axis.
%           The very first and the very last bin are special: they are the
%           sticky decision bounds.
%
% frbins    Vector indexing values of normalized firing rates
%
% Pjoints   Matrix numel(t0s)-x-numel(frbins)-x-numel(x) in size. Each of
%           the t0 entries corresponds to the joint probability
%           distribution P(r, a | t=t0).
%
% fr_given_as   Marix numel(t0s)-x-numel(x) in size. Each row is a vector
%          for the corresponding t0 value, and represents the mean firing
%          rate given the value of a.
%
% fr_var_given_as Marix numel(t0s)-x-numel(x) in size. Each row is a vector
%          for the corresponding t0 value, and represents the variance of the firing
%          rate given the value of a.
%
% a_given_frs   Marix numel(t0s)-x-numel(frbins) in size. Each row is a vector
%          for the corresponding t0 value, and represents the mean value of
%          given a firing rate. 
%
% To display the joint probability distribution (for example, for the 2nd
% timepoint t0), we recommend 
%
%   imagesc(x(2:end-1), frbins(2:end), Pjoints(2, 2:end, 2:end-1));
%   colormap hot; xlabel('a'); ylabel('norm. firing rate');
%
% The reason for cutting off the first timepoint is that the zero firing
% rate bin often has a lot of mass, and will thus be off scale w.r.t. color axis; 
% same for cutting off the sticky bound bins, since they also accumulate
% more probability mass than typical bins and therefore throw the scale
% off.
%

function [constant_x, frbins, Pjoints, fr_given_as, fr_var_given_as, a_given_frs] = dyn_compile_dv2(cellid, t0s, p, varargin)

frbins=0;
pairs = { ...
	'lag'        0.2         ; ...
	'frbins'     0:0.1:50    ; ...
	'dt'         0.01        ; ...
    'alignment'      'stimstart-cout-mask' ; ...   % string matching one of align_strs in align_LUT
	'direction'       'backward' ; ...   % 'forward' or 'backward'
	'trialnums'  []          ; ...
	'use_nans'   0           ; ...
    'krn_width'  []          ; ...      % forces neural data to have this kernel width; if empty, uses whatever is in data    'krn_type'   []          ; ...      % forces neural data to have this kernel type; if empty, uses whatever is in data
    'krn_type'   'halfgauss' ; ...      % forces neural data to have this kernel type
    'fr_dt'      []          ; ...      % forces neural data to have this bin size; if empty, uses whatever is in data
    'norm_type'     'div'    ; ...      % type of firing rate normalization; 'div' is divisive, 'z' is z-score, 'none' no normalization
    'fit_version'   'byrat'   ; ...      % fit params to use; 'byrat', 'combined'
    'n_iter'    1            ; ...      % number of iterations for refining estimate of DV
    'param_scale_num'        1     ; ... % parameter number to scale
    'param_scale_factor'     1     ; ... % multiplicative factor of that parameter
    'average_time_points'   0;
}; parseargs(varargin, pairs);

% check for bad inputs on fit_version
if ~strcmp(fit_version, 'byrat')
    error('This code isnt set up for meta rat models')
end

% Load data for this rat
if cellid < 0
    % This is a fake neuron for testing purposes;
    disp('Synthetic Neuron')
    data = dyn_cell_packager('17784'); % 17784/H066/565685, 16821/H037, 18548/H129, 18417/H084, 18532/H066/601406
else
    data = dyn_cell_packager(cellid);
end


% make sure krn_width, krn_type and fr_dt are correct; otherwise re-package with correct values
repack = false;
if ~isempty(krn_width) && data.krn_width~=krn_width
    repack = true;
else
    krn_width = data.krn_width;
end
if ~isempty(fr_dt) && data.bin_size~=fr_dt
    repack = true;
else
    fr_dt = data.bin_size;
end
if ~isfield(data,'krn_type')    % remove this part for after re-packaging all data
    repack = true;
elseif ~isempty(krn_type) && ~strcmp(data.krn_type,krn_type)
    repack = true;
else
    krn_type = data.krn_type;
end
if repack
    if cellid < 0
        error('repacking for a synthetic neuron...that shouldnt happen');
    end
    data = dyn_cell_packager(cellid,'krn_width',krn_width,'krn_type',krn_type,'bin_size',fr_dt,'force',1);
end

% Check to make sure we have the right number of trials
if ~isempty(trialnums),
	trials = trialnums;
else
	trials = 1:data.ngood_trials;
end;

% Set up some storage
Pjoints         = [];
fr_given_as     = [];
fr_var_given_as = [];
a_given_frs     = [];
total_mass      = zeros(numel(t0s),1);

% find alignment index
[align_strs, align_args] = align_LUT;
align_ind = strmatch(alignment,align_strs,'exact');

% find the offset between model time and fr time based on the alignment
% model is always aligned in time to stim_start; fr is not
% get reference event
ref_event = get_ref_event(align_args{align_ind});
try
    ref_times = eval(['data.trials.' ref_event]);
catch
    error('Error while evaluating ref_event string.')
end
ss_times = data.trials.stim_start; % stim_start is ref event for model
mt_offset = ref_times - ss_times;  % offset between model and fr_t on each trial

% load model posterior for all trials
load(['/home/alex/Dropbox/spikes/model/model_posterior_' num2str(data.sessid) '.mat']);

% Need to filter trials from model to match data variable
p.reload    = 0;
p.ratname   = data.ratname;
[~, vec_data, ~,~] = get_behavior_data('/home/alex/Dropbox/spikes/data/', data.cellid, data.sessid,p);
model = model(vec_data.good);


% do some sanity checks to make sure I aligned the model responses correctly. 
if numel(trials) ~= numel(model)
    error('Wrong number of trials in model, are you clearing bad trials consistently?')
elseif sum(data.trials.gamma == vec_data.gamma) ~= numel(trials)
    error('Bad trial alignment between model predictions and cell packager');
end

% Check for variable a-binning, which happens with large sensory noise to save space
numabins = length(model(1).posterior.avals);
constant_x = model(1).posterior.avals;
allsame=true;
for ti=2:numel(trials)
    if length(model(ti).posterior.avals) ~= numabins
        allsame=false;
        if length(model(ti).posterior.avals) > numabins
            numabins=length(model(ti).posterior.avals);
            constant_x = -10:0.1:10;
        end
    end
end

% Iterate over trials
if (cellid > -11 & cellid <= -6) | cellid == -4
    vars = zeros(1,numel(trials));
end
for ti=1:numel(trials)
    trialnum = trials(ti);
    
    % get this trials model predictions
    x = model(ti).posterior.avals;
    t = model(ti).posterior.T;
    P = model(ti).posterior.pdf;

    % Figure out how to normalize firing rate
    ft = data.frate_t{align_ind};
    switch norm_type
        case 'none'
            fr = data.frate{align_ind}(trialnum,:);
        case 'div'
            fr = data.frate{align_ind}(trialnum,:) ./ data.norm_mean;
        case 'z'
            fr = (data.frate{align_ind}(trialnum,:) - data.norm_mean) ./ data.norm_std;
        otherwise
            error('Do not recognize norm_type')
    end
    if cellid < 0
        % this is a fake neuron for testing purposes, replace actual fr with synthetic fr
        fr = get_synthetic_fr(cellid, fr,ft,model(ti).posterior.mean,norm_type);

        if cellid == -6
            for jj = 1:size(P,1)
                ps = sum(P(jj,:));
                temp = P(jj,:);
                temp = temp.^2;
                temp = temp./sum(temp);
                temp = temp.*ps;
                P(jj,:) = temp;
            end
        elseif cellid == -7
            for jj = 1:size(P,1)
                ps = sum(P(jj,:));
                temp = P(jj,:);
                temp = temp.^4;
                temp = temp./sum(temp);
                temp = temp.*ps;
                P(jj,:) = temp;
            end
        elseif cellid == -8
            for jj = 1:size(P,1)
                ps = sum(P(jj,:));
                temp = P(jj,:);
                temp = temp.^6;
                temp = temp./sum(temp);
                temp = temp.*ps;
                P(jj,:) = temp;
            end
        elseif cellid == -9
            for jj = 1:size(P,1)
                ps = sum(P(jj,:));
                temp = P(jj,:);
                temp = temp.^8;
                temp = temp./sum(temp);
                temp = temp.*ps;
                P(jj,:) = temp;
            end
        elseif cellid == -10
            for jj = 1:size(P,1)
                ps = sum(P(jj,:));
                temp = P(jj,:);
                temp = temp.^10;
                temp = temp./sum(temp);
                temp = temp.*ps;
                P(jj,:) = temp;
            end
        end
        if cellid == -4 | (cellid <= -6 & cellid > -11)

            trial_vars = zeros(1,size(P,1));
            for jj=1:size(P,1);
                temp = P(jj,:);
                temp = temp./sum(temp);
                thismean = sum(temp.*x);
                varx = (x - thismean).^2;
                trial_vars(jj) = sum(temp.*varx); 
            end
            vars(ti) = mean(trial_vars);
        end
    end

    
    % Check if time bins on firing rate is too coarse relative to tuning curve time bins
    if min(diff(t0s)) < min(diff(ft))
        error('Time bins for tuning curves are smaller than firing rate smoothed vector. Recompute firing rate.')
    end

    % Check to see if the time bins requested for the tuning curve extend beyond the firing rate vector
    if t0s(end) > ft(end)
        error('The last time point is beyond the firing rate time points. Decrease the time bins t0s, or recompute the firing rate')
    end
    
 
    % Iterate over time indexes
    my_fr = [];
    for time_i=1:numel(t0s),
        t0 = t0s(time_i);
       
        % Set up storage 
        if isempty(Pjoints),
            if ~allsame
                numx = 201;
            else
                numx = numel(x);
            end
            if use_nans,
                Pjoints     = NaN*ones(numel(t0s), numel(frbins), numx);
            else
                Pjoints     = zeros(numel(t0s), numel(frbins), numx);
            end;
            fr_given_as     = zeros(numel(t0s), numx);
            fr_var_given_as = zeros(numel(t0s), numx);
            a_given_frs     = zeros(numel(t0s), numel(frbins));
        end
     
        mask_after_stim_firing  = p.mask_after_stim_firing; 
        mask_stim_firing        = p.mask_stim_firing; 
        average_a_vals          = p.average_a_vals;
        average_fr              = p.average_fr;

        if mask_stim_firing && mask_after_stim_firing
            error('You cant mask stimulus firing and after-stimulus firing')
        end

        if sum(mt_offset ~= 0) > 0 
            error('Be Careful, Im not sure mt_offset was implemented correctly with time averaging')
        end
        if lag > 0 
            error('I havent implemented lag with time averaging yet')
        end
    
        % Compute which time bins to average together. 
        % Take a bin centered at t0, and look at the half-interval between the previous and next t0.
        t0_durs     = diff(t0s);
        if time_i == numel(t0s)
            this_dur    = 0;
        else
            this_dur    = t0_durs(time_i);
        end
        if time_i > 1
            last_dur = t0_durs(time_i - 1);
        else
            last_dur = 0;
        end
        if (this_dur/2 + t0) > ft(find(~isnan(fr),1,'last'))
            this_dur = 2*(ft(find(~isnan(fr),1,'last')) - t0);
        end
        all_time_indexes = t > (t0-last_dur/2) & t < (t0 + this_dur/2);
        
        % Average over a-vals from multiple time bins, or select just the current bin
        if average_a_vals
            time_indexes = all_time_indexes;
        else
            [~, time_indexes]   = min(abs(t-t0)); 
        end
    
        % Average over firing rates from multiple time bins, or select just the current bin
        r0_mesh = interp1(ft, fr, [t0-last_dur/2, t0, t0+this_dur/2]);
        if average_fr
            r0 = mean(r0_mesh);% The mean firing rate during the time window
        else
            r0 = interp1(ft, fr, t0+lag); 
        end
        [~, fbin] = min(abs(frbins-r0));     % fbin is our bin in frbins

        % Do you want to mask spiking activity after the stimulus stopped?
        % Should you include this firing rate? If its NaN, or if we are masking post-stimulus activity
        if mask_after_stim_firing
            include_fr = ~isnan(r0) & (t(end)-t0-mt_offset(trialnum) > 0); 
        elseif mask_stim_firing
            include_fr = ~isnan(r0) & (t(end)-t0-mt_offset(trialnum) < 0); 
        else
            include_fr = ~isnan(r0);
        end
       
        % Add this bin to Pjoints
        if include_fr
            if ~allsame
                thisx = size(P,2); 
                sd = ceil(thisx/2) - 100;
                ed = ceil(thisx/2) + 100;
                this_prob = mean(P(time_indexes,sd:ed),1);
            else
                this_prob = mean(P(time_indexes,:),1);
            end
            this_prob = (this_prob./sum(this_prob)).*10;
            if use_nans && isnan(Pjoints(time_i, fbin, 2))  % this should only happen is use_nans was 1
                Pjoints(time_i, fbin, :)    = this_prob;%mean(P(time_indexes,sd:ed),1);
            else
                Pjoints(time_i, fbin, :)    = squeeze(Pjoints(time_i, fbin,:))' + this_prob;%mean(P(time_indexes,:),1); %#ok<*AGROW>
            end
            total_mass(time_i)              = total_mass(time_i) + 1;
        end

        if 0 % DEBUG CODE for averaging firing rates
            clf; hold on;
            plot(ft,fr, 'bo-')
            plot(ft(ft_indexes), fr(ft_indexes), 'ro','markerfacecolor','r')
            plot([t0-last_dur/2, t0, t0+this_dur/2], r0_mesh, 'mx','markersize',15)
            plot([t0-last_dur/2, t0+this_dur/2], [r0 r0], 'm-')
            plot(t0s(1:time_i), my_fr, 'k-')
            keyboard;
        end
        if 0 % DEBUG CODE for averaging A-vals
            clf; hold on;
            temp = P(time_indexes,:);
            temp1 = mean(temp,1);
            plot(x,temp)
            plot(x, temp1, 'k-','linewidth',4)
            keyboard;
        end
    end;
end;
if (cellid <= -6 & cellid > -11) | cellid == -4
    %disp(mean(vars))
    average_var = mean(vars);
    save(['deconvolve/var_data_' num2str(abs(cellid)) '.mat'],'average_var'); 
end


% Not sure exactly what is happening here 
for time_i=1:numel(t0s)
    Pjoints(time_i,:,:) = Pjoints(time_i,:,:)/total_mass(time_i);
    
    myPj = squeeze(Pjoints(time_i,:,:));
    
    Pjoint_given_a = myPj ./ (ones(size(myPj,1),1)*sum(myPj,1));
    Pjoint_given_f = myPj ./ (sum(myPj,2)*ones(1,size(myPj,2)));
    
    fr_given_as(time_i,:)       = Pjoint_given_a'*frbins';
    fr_var_given_as(time_i,:)   = (Pjoint_given_a'*(frbins'.^2)) - (fr_given_as(time_i,:)'.^2);
    a_given_frs(time_i,:)       = Pjoint_given_f*constant_x';
end;

% get a string of the ref_event in the arguments that come from the chosen
% alignment
function ref_event = get_ref_event(arguments)
ind = 1;  
ref_event = 'cpoke_end';    % default alignment in make_rate_functions
while ind <= length(arguments),
    if strmatch(arguments{ind},'ref_event','exact')
        ref_event = arguments{ind+1};
        return
    end 
    ind = arg+2; 
end;
