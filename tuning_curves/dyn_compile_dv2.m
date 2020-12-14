function [constant_x, frbins, Pjoints, fr_given_as, fr_var_given_as, a_given_frs]...
    = dyn_compile_dv2(cellid, t0s, ops, varargin)

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
% alignment string matching one of align_strs in dyn_align_LUT. This will set
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



frbins=0;
p = inputParser;
addParameter(p,'lag',        0.2         );
addParameter(p,'frbins',     0:0.1:50    );
addParameter(p,'dt',         0.01        );
addParameter(p,'alignment',      'stimstart-cout-mask' );   % string matching one of align_strs in dyn_align_LUT
addParameter(p,'direction',       'backward' );   % 'forward' or 'backward'
addParameter(p,'trialnums',  []          );
addParameter(p,'use_nans',   0           );
addParameter(p,'krn_width',  []          );      % forces neural data to have this kernel width; if empty, uses whatever is in data    'krn_type'   []          );      % forces neural data to have this kernel type; if empty, uses whatever is in data
addParameter(p,'krn_type',   'halfgauss' );      % forces neural data to have this kernel type
addParameter(p,'fr_dt',      []          );      % forces neural data to have this bin size; if empty, uses whatever is in data
addParameter(p,'norm_type',     'div'    );      % type of firing rate normalization; 'div' is divisive, 'z' is z-score, 'none' no normalization
addParameter(p,'fit_version',   'byrat'   );      % fit params to use; 'byrat', 'combined'
addParameter(p,'n_iter',    1            );      % number of iterations for refining estimate of DV
addParameter(p,'param_scale_num',        1     ); % parameter number to scale
addParameter(p,'param_scale_factor',     1     ); % multiplicative factor of that parameter
addParameter(p,'average_time_points',   0);
parse(p,varargin{:})
struct2vars(p.Results);

dp = set_dyn_path;
spike_data_dir = dp.spikes_dir;

% find alignment index
[align_strs, align_args] = dyn_align_LUT;
align_ind = strmatch(alignment,align_strs,'exact');

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
[repack, krn_width, fr_dt, krn_type] = check_params(data,krn_width,fr_dt,krn_type); %#ok<NODEF>
if repack
    if cellid < 0
        error('repacking for a synthetic neuron...that shouldnt happen');
    end
    data = dyn_cell_packager(cellid,'krn_width',krn_width,'krn_type',krn_type,'bin_size',fr_dt,'force',1);
end

% get some info about the trials
ops.reload    = 0;
ops.ratname   = data.ratname;
[~, vec_data, ~,~] = get_behavior_data(spike_data_dir, data.cellid, data.sessid,ops);

[model, constant_x, allsame] = get_data_model_p(data, vec_data);

% Check to make sure we have the right number of trials
if ~isempty(trialnums)
    trials = trialnums;
else
    trials = 1:data.ngood_trials;
end;

% do some sanity checks to make sure I aligned the model responses correctly.
if numel(trials) ~= numel(model)
    error('Wrong number of trials in model, are you clearing bad trials consistently?')
elseif sum(data.trials.gamma == vec_data.gamma) ~= numel(trials)
    error('Bad trial alignment between model predictions and cell packager');
end

% Iterate over trials
if (cellid > -11 & cellid <= -6) | cellid == -4
    vars = zeros(1,numel(trials));
end

% find the offset between model time and fr time based on the alignment
% model is always aligned in time to stim_start; fr is not
% get reference event
ref_ev_ind  = find(cellfun(@(x) strcmp('ref_event',x),align_args{align_ind})) + 1;
ref_event   = align_args{align_ind}{ref_ev_ind};
ref_times   = data.trials.(ref_event);
ss_times    = data.trials.stim_start; % stim_start is ref event for model
mt_offset   = ref_times - ss_times;  % offset between model and fr_t on each trial

if sum(mt_offset ~= 0) > 0
    error('Be Careful, Im not sure mt_offset was implemented correctly with time averaging')
end

% Set up some storage
Pmass           = [];
total_mass      = zeros(numel(t0s),1);

for ti=1:numel(trials)
    trialnum = trials(ti);

    % get this trials model predictions
    x = model(ti).posterior.avals;
    t = model(ti).posterior.T;
    P = model(ti).posterior.pdf;
    
    if isempty(Pmass)
        if ~allsame
            numx = 201;
        else
            numx = numel(x);
        end
        if use_nans
            Pmass     = NaN*ones(numel(t0s), numel(frbins), numx);
        else
            Pmass     = zeros(numel(t0s), numel(frbins), numx);
        end;
    end
    
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
        [fr, P] = synthetic_fr_P(cellid, fr, ft, model, norm_type)
    end
    
    % Check if time bins on firing rate is too coarse relative to tuning curve time bins
    if min(diff(t0s)) < min(diff(ft))
        error('Time bins for tuning curves are smaller than firing rate smoothed vector. Recompute firing rate.')
    end
    
    % Check to see if the time bins requested for the tuning curve extend beyond the firing rate vector
    if t0s(end)+lag > ft(end)
        error('The last time point is beyond the firing rate time points. Decrease the time bins t0s, or recompute the firing rate')
    end
    
    % Iterate over time indexes
    mask_after_stim_firing  = ops.mask_after_stim_firing;
    mask_stim_firing        = ops.mask_stim_firing;
    average_a_vals          = ops.average_a_vals;
    average_fr              = ops.average_fr;
    if mask_stim_firing && mask_after_stim_firing
        error('You cant mask stimulus firing and after-stimulus firing')
    end
    my_fr       = [];
    t0_durs     = [diff(t0s) 0];

    for time_i=1:numel(t0s)
        t0 = t0s(time_i);
        
        % Compute which time bins to average together for A values
        % Take a bin centered at t0, and look at the half-interval between the previous and next t0.
        if time_i == 1
            last_dur = 0;
        else
            last_dur = t0_durs(time_i - 1);
        end
        this_dur    = t0_durs(time_i);
        
        if (this_dur/2) > ft(find(~isnan(fr),1,'last'))
            this_dur_nolag = 2*(ft(find(~isnan(fr),1,'last')) - t0);
        else
            this_dur_nolag = this_dur;
        end
        if (this_dur/2 + t0+lag) > ft(find(~isnan(fr),1,'last'))
            this_dur_lag = 2*(ft(find(~isnan(fr),1,'last')) - t0-lag);
        else
            this_dur_lag = this_dur;
        end
        
        all_time_indexes = t > (t0-last_dur/2) & t < (t0 + this_dur_nolag/2);
        
        % Average over a-vals from multiple time bins, or select just the current bin
        if average_a_vals
            time_indexes = all_time_indexes;
        else
            [~, time_indexes]   = min(abs(t-t0));
        end

        % Average over firing rates from multiple time bins, or select just the current bin
        r0_mesh = interp1(ft, fr, [t0+lag-last_dur/2, t0+lag, t0+lag+this_dur_lag/2]);
        if average_fr
            r0 = mean(r0_mesh);% The mean firing rate during the time window
        else
            r0 = interp1(ft, fr, t0+lag);
        end
        [~, fbin] = min(abs(frbins-r0));     % fbin is our bin in frbins
        
        % Do you want to mask spiking activity after the stimulus stopped?
        % Should you include this firing rate? If its NaN, or if we are masking post-stimulus activity
        if mask_after_stim_firing
            include_fr = ~isnan(r0) & (t(end)-(t0+lag)-mt_offset(trialnum) > 0);
        elseif mask_stim_firing
            include_fr = ~isnan(r0) & (t(end)-(t0+lag)-mt_offset(trialnum) < 0);
        else
            include_fr = ~isnan(r0);
        end
        
        % Add this bin to Pmass
        if include_fr
            if ~allsame
                thisx = size(P,2);
                sd = ceil(thisx/2) - 100;
                ed = ceil(thisx/2) + 100;
                this_prob = mean(P(time_indexes,sd:ed),1);
            else
                this_prob = mean(P(time_indexes,:),1);
            end
            this_prob = (this_prob./sum(this_prob)).*1;
            if use_nans && isnan(Pmass(time_i, fbin, 2))  % this should only happen is use_nans was 1
                Pmass(time_i, fbin, :)    = this_prob; %#ok<*AGROW>
                %mean(P(time_indexes,sd:ed),1);
            else
                Pmass(time_i, fbin, :)    = squeeze(Pmass(time_i, fbin,:))' + this_prob; %#ok<*AGROW>
                %mean(P(time_indexes,:),1); 
            end
            total_mass(time_i)              = total_mass(time_i) + 1; % keep track of how many trials contribute to each time bin
        end
        if 0 &time_i == numel(t0s) % DEBUG CODE for averaging firing rates
            clf; hold on;
            plot(ft,fr, 'bo-')
            %     plot(ft(ft_indexes), fr(ft_indexes), 'ro','markerfacecolor','r')
            plot([t0+lag-last_dur/2, t0+lag, t0+lag+this_dur/2], r0_mesh, 'mx','markersize',15)
            plot([t0+lag-last_dur/2, t0+lag+this_dur/2], [r0 r0], 'm-')
            my_fr = [my_fr r0];
            plot(t0s(1:time_i)+lag, my_fr, 'k-')
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


% get joint and conditional distributions 
fr_given_as     = zeros(numel(t0s), numx);
fr_var_given_as = zeros(numel(t0s), numx);
a_given_frs     = zeros(numel(t0s), numel(frbins));
Pjoints = nan(size(Pmass));
for time_i=1:numel(t0s)
    Pjoints(time_i,:,:) = Pmass(time_i,:,:)/total_mass(time_i); % so far Pjoints is a sum, so we're turning into a probability here
    
    myPj = squeeze(Pjoints(time_i,:,:)); 
    
    Pjoint_given_a = myPj ./ (ones(size(myPj,1),1)*sum(myPj,1)); % normalize Pjoint to get firing rate tuning curve wrt a
    Pjoint_given_f = myPj ./ (sum(myPj,2)*ones(1,size(myPj,2))); % normalize Pjoint to get a distribution wrt firing rate 
    
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



function [fr, P] = synthetic_fr_P(cellid, fr, ft, model, norm_type)
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

function [repack, krn_width, fr_dt, krn_type] = check_params(data,krn_width,fr_dt,krn_type)
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
