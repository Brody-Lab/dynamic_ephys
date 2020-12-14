function [constant_x, frbins, Pjoints, fr_given_as, fr_var_given_as, a_given_frs] = ...
    fr_dv_switch(cellid, t0s, ops, varargin)
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
addParameter(p,'which_switch',   'model');
addParameter(p,'array_data',   []);
addParameter(p,'vec_data',   []);
addParameter(p,'clear_bad_strengths',   1);
addParameter(p,'bad_strength',   0);
addParameter(p,'fit_line',   1);
addParameter(p,'final_only',   0);
addParameter(p,'exclude_final',   0);

parse(p,varargin{:});
struct2vars(p.Results);

if nargin < 1
    cellid = 18181;
end

dp = set_dyn_path;
spike_data_dir = dp.spikes_dir;

data            = dyn_cell_packager(cellid);

% load and clean up array_data and vec_data for this session
% load and clean up array_data and vec_data for this session
[switch_to_0, switch_to_1, array_data,vec_data] = ...
    get_switches(cellid, 'which_switch',which_switch,...
    'array_data',array_data,'vec_data',vec_data,...
    'clear_bad_strengths', clear_bad_strengths, ...
    'bad_strength', bad_strength, 'fit_line', fit_line,...
    'exclude_final', exclude_final, 'final_only', final_only);

% decide which trials to analyze and which to discard
has_switch  = (cellfun(@length, switch_to_0) + cellfun(@length, switch_to_1)) > 0;
includes    = has_switch(:);
n_trials    = sum(includes);

%switch_to_0 = switch_to_0(includes);
%switch_to_1 = switch_to_1(includes);

n_switch_to_0 = cellfun(@length, switch_to_0);
n_switch_to_1 = cellfun(@length, switch_to_1);

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
[align_strs, align_args] = dyn_align_LUT;
align_ind   = strmatch(alignment,align_strs,'exact');
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
counter         = 0;
for ti=1:numel(trials)
    trialnum = trials(ti);

    % get this trials model predictions
    x = model(ti).posterior.avals;
    t = model(ti).posterior.T;
    P = model(ti).posterior.pdf;
    
    this_switches   = sort([switch_to_0{ti} switch_to_1{ti}]);
    this_nswitches  = length(this_switches);   
    
    % Iterate over time indexes
    mask_after_stim_firing  = ops.mask_after_stim_firing;
    mask_stim_firing        = ops.mask_stim_firing;
    average_a_vals          = ops.average_a_vals;
    average_fr              = ops.average_fr;
    if mask_stim_firing && mask_after_stim_firing
        error('You cant mask stimulus firing and after-stimulus firing')
    end
    
    for si = 1:this_nswitches
        counter     = counter + 1;
        switch_t0s  = this_switches(si) + t0s; 
        
        if si == this_nswitches
            mask_after  = array_data(ti).stim_end;
        else
            mask_after  = this_switches(si+1);
        end
        
        if si == 1
            mask_before  = 0;
        else
            mask_before  = this_switches(si-1);
        end
        good_tinds = find(switch_t0s > mask_before & switch_t0s < mask_after);
        %switch_t0s = switch_t0s(good_tinds);
        
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
        if (min(diff(ft) - min(diff(switch_t0s)))) > .001
            error(['Time bins for tuning curves are smaller than' ...
                'firing rate smoothed vector. Recompute firing rate.']);
        end
        
        % Check to see if the time bins requested for the tuning curve extend beyond the firing rate vector
        if switch_t0s(good_tinds(end))+lag > (ft(end) + min(diff(ft)))
            error(['The last time point is beyond the firing rate' ...
                'time points. Decrease the time bins t0s, or recompute' ...
                'the firing rate']);
        end
        

        my_fr       = [];
        t0_durs     = [diff(switch_t0s) 0];
               
        for time_i = good_tinds
            t0 = switch_t0s(time_i);
            
            % Compute which time bins to average together for A values
            % Take a bin centered at t0, and look at the half-interval between the previous and next t0.
            if time_i == good_tinds(1)
                last_dur = 0;
            else
                last_dur = t0_durs(time_i - 1);
            end
            if time_i == good_tinds(end)
                this_dur = 0;
            else
                this_dur = t0_durs(time_i);
            end
            
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
            if 0 &time_i == numel(switch_t0s) % DEBUG CODE for averaging firing rates
                clf; hold on;
                plot(ft,fr, 'bo-')
                %     plot(ft(ft_indexes), fr(ft_indexes), 'ro','markerfacecolor','r')
                plot([t0+lag-last_dur/2, t0+lag, t0+lag+this_dur/2], r0_mesh, 'mx','markersize',15)
                plot([t0+lag-last_dur/2, t0+lag+this_dur/2], [r0 r0], 'm-')
                my_fr = [my_fr r0];
                plot(switch_t0s(1:time_i)+lag, my_fr, 'k-')
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
        
    end
    


    
end;
if (cellid <= -6 & cellid > -11) | cellid == -4
    %disp(mean(vars))
    average_var = mean(vars);
    save(['deconvolve/var_data_' num2str(abs(cellid)) '.mat'],'average_var');
end


% get joint and conditional distributions 
fr_given_as     = zeros(numel(switch_t0s), numx);
fr_var_given_as = zeros(numel(switch_t0s), numx);
a_given_frs     = zeros(numel(switch_t0s), numel(frbins));
Pjoints = nan(size(Pmass));
for time_i=1:numel(switch_t0s)
    Pjoints(time_i,:,:) = Pmass(time_i,:,:)/total_mass(time_i); % so far Pjoints is a sum, so we're turning into a probability here
    
    myPj = squeeze(Pjoints(time_i,:,:)); 
    
    Pjoint_given_a = myPj ./ (ones(size(myPj,1),1)*sum(myPj,1)); % normalize Pjoint to get firing rate tuning curve wrt a
    Pjoint_given_f = myPj ./ (sum(myPj,2)*ones(1,size(myPj,2))); % normalize Pjoint to get a distribution wrt firing rate 
    
    fr_given_as(time_i,:)       = Pjoint_given_a'*frbins';
    fr_var_given_as(time_i,:)   = (Pjoint_given_a'*(frbins'.^2)) - (fr_given_as(time_i,:)'.^2);
    a_given_frs(time_i,:)       = Pjoint_given_f*constant_x';
    
end;


function fr_dv_switch_(fr, switch_ts)

