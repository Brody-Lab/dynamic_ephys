function res = compute_switch_triggered_average(data,varargin)
% function res =  compute_switch_triggered_average(data, varargin)
%
% Computes a switch triggered average for cells recording in the dynamic
% clicks task.
%
% input:
% -----
% data        - a struct produced by calling package_dyn_phys for a given
%               cellid. if input is numeric, code assumes it is a cellid
%               and passes it to package_dyn_phys
%
% output:
% -----
% res a structure containing switch triggered responses and d' for the real
% and, if desired, shuffled data
%
% example usage:
% res = compute_switch_triggered_average(17784, 'n_shuffles', 1000,
%                       'array_data', adata, 'vec_data', vdata)
% computes sta for cell 17784 using the model switches contained in
% adata.model_switch_to_0, also recomputes sta 1000 times using
% shuffles of the switch times across trials
%
%
% res = compute_switch_triggered_average(data, 'clear_bad_strengths', 0)
% Computes sta using the data struct. It computes the model switch times,
% but does not remove low salience switches, removing only the first model
% switch (which is really a switch from agnosticism). No shuffles are
% performed.

% optional parameters specified by user
p = inputParser;
addParameter(p, 'n_shuffles', 0) % whether to run a shuffle test
addParameter(p, 'which_switch', 'model') % other option is 'generative'
addParameter(p, 'align_ind', 1) % 1 means align to stim on and masked after stim off
addParameter(p, 'norm_type', 'none') % 'div', 'log_div', 'z' are other options
addParameter(p, 'post', 2) % time over which to make this plot
addParameter(p, 'model_smooth_wdw', 100);
addParameter(p, 'array_data', []);
addParameter(p, 'vec_data', []);
addParameter(p, 'include_str', 'data.trials.hit == 1');
addParameter(p, 'clear_bad_strengths', 1);
addParameter(p, 'bad_strength', 0);
addParameter(p, 'save_dir', '');
addParameter(p, 'model_dir', '');
addParameter(p, 'save_file', 1);
addParameter(p, 'force', 0); % determines whether to recompute
addParameter(p, 'mask_other_switch', 0); % determines whether to recompute
addParameter(p, 'fit_line', 1); % determines method for quantifying state strength
addParameter(p, 'condition_residual', 0); 
addParameter(p, 'compute_residual', 1);
addParameter(p, 'exclude_final', 0); 
addParameter(p, 'final_only', 0); 
parse(p,varargin{:});
p = p.Results;

dp = set_dyn_path;
if isempty(p.save_dir)
    p.save_dir = dp.celldat_dir;
end
if isempty(p.model_dir)
    p.model_dir = dp.model_mean_dir;
end

% decide how to interpret first input
if isnumeric(data)
    cellid = data;
else
    cellid = data.cellid;
end


if (~isempty(p.save_dir) && ~exist(p.save_dir, 'dir')) || ...
        (~isempty(p.model_dir) && ~exist(p.model_dir, 'dir'))
    error('your default save or model directory don''t exist')
end

% determine the file path for this STA and load the file if desired
if ~isempty(p.save_dir)
    if p.clear_bad_strengths
        save_path = fullfile(p.save_dir, [p.which_switch '_STA_' num2str(cellid) '_' num2str(p.bad_strength) '.mat']);
    else
        save_path = fullfile(p.save_dir, [p.which_switch '_STA_' num2str(cellid) '.mat']);
    end
    if ~p.force & exist(save_path,'file')
        load(save_path,'res');
        return
    elseif ~exist(p.save_dir,'dir') | isempty(p.save_dir)
        warning('save dir incorrectly specified')
    end
end

disp('computing switch triggered average')
% if data is a cellid, load it
if isnumeric(data)
    data = dyn_cell_packager(data);
end

% load and clean up array_data and vec_data for this session
[switch_to_0, switch_to_1, p.array_data, p.vec_data] = ...
    get_switches(cellid, 'which_switch',p.which_switch,...
    'array_data',p.array_data,'vec_data',p.vec_data,...
    'clear_bad_strengths', p.clear_bad_strengths, ...
    'bad_strength', p.bad_strength, 'fit_line', p.fit_line,...
    'exclude_final', p.exclude_final, 'final_only', p.final_only);

% decide which trials to analyze and which to discard
has_switch  = (cellfun(@length, switch_to_0) + cellfun(@length, switch_to_1)) > 0;
includes    = eval(p.include_str) & has_switch';
n_trials    = sum(includes);

switch_to_0 = switch_to_0(includes);
switch_to_1 = switch_to_1(includes);

n_switch_to_0 = cellfun(@length, switch_to_0);
n_switch_to_1 = cellfun(@length, switch_to_1);

% normalize firing rates
switch p.norm_type
    case 'none'
        norm_ys = data.frate{p.align_ind};
    case 'div'
        norm_ys = data.frate{p.align_ind} ./ data.norm_mean;
    case 'log_div'
        norm_ys = log(data.frate{p.align_ind} ./ data.norm_mean);
    case 'z'
        norm_ys = (data.frate{p.align_ind} - data.norm_mean) ./ data.norm_std;
    otherwise
        error('norm_type not recognized')
end

% time axis
t               = data.frate_t{p.align_ind};
[~, pre_bin]    = min(abs(t));
[~, post_bin]   = min(abs(t-p.post));
t               = t(pre_bin:post_bin);
% lags are time differences between switches and spikes
lags    = [-fliplr(t(pre_bin+1:post_bin)) t(pre_bin:post_bin)];
n_lags  = length(lags);

% initialize switch triggered response for matrix for each side
STR_right = nan(sum(n_switch_to_1), n_lags, p.n_shuffles+1);
STR_left = nan(sum(n_switch_to_0), n_lags, p.n_shuffles+1);

% get data from included trials and relevant time bins
norm_ys = norm_ys(includes,pre_bin:post_bin);
right_trials = data.trials.rat_dir(includes)==1;
left_trials = ~right_trials;

% calculate mean firing rate
fr_mean = nanmean(norm_ys);

% compute residual firing rates by subtracting the mean from all trials
if ~p.compute_residual
    fr_residual = norm_ys;
elseif p.condition_residual
    fr_residual = norm_ys;
    left_mean   = nanmean(norm_ys(left_trials,:));
    right_mean  = nanmean(norm_ys(right_trials,:));
    fr_residual(left_trials,:) = norm_ys(left_trials,:) - left_mean;
    fr_residual(right_trials,:) = norm_ys(right_trials,:) - right_mean
    if p.n_shuffles
        error('not sure if we are handling shuffle test with a choice-conditioned residual')
    end
else
    fr_residual = norm_ys - repmat(fr_mean,n_trials,1);
end
% if desired, shuffle residual firing rates across included trials
real_ind = 1:size(norm_ys,1);
test_ind = repmat(real_ind,p.n_shuffles+1,1);
for ss = 2:p.n_shuffles
    test_ind(ss,:) = real_ind(randperm(length(real_ind)));
end

% loop through trials to build the sta matrix
n_computed_switches_to_0 = 0;
n_computed_switches_to_1 = 0;
for tt=1:n_trials
    % make response vector
    s2_full = fr_residual(test_ind(:,tt),:)';
    
    this_switch_1 = switch_to_1{tt};
    this_switch_0 = switch_to_0{tt};
    this_switch = [this_switch_1 this_switch_0];
    this_state  = [ones(size(this_switch_1)) zeros(size(this_switch_0))];
    [this_switch, sort_dex] = sort(this_switch);
    this_state = this_state(sort_dex);
    s2 = s2_full;

    if p.mask_other_switch
        for m = 1:length(this_switch)
            [t_start, t_end, ind] = get_lag_params(t,this_switch(m),n_lags);
            s2 = s2_full;
            if m > 1
                s2(t<this_switch(m-1),:) = nan;
            end
            if m < length(this_switch)
                s2(t>this_switch(m+1),:) = nan;
            end
            
            if this_state(m) == 0
                n_computed_switches_to_0 = n_computed_switches_to_0+1;
                if sum(ind)>0
                    STR_left(n_computed_switches_to_0,max(1,t_start):t_end,:) = s2(ind,:);
                end
            elseif this_state(m) == 1
                n_computed_switches_to_1 = n_computed_switches_to_1+1;
                if sum(ind)>0
                    STR_right(n_computed_switches_to_1,max(1,t_start):t_end,:) = s2(ind,:);
                end
            end
        end
    else
        
        %loop over switches to state 0 and align this spike train to the right time lags
        if ~isempty(switch_to_0{tt})
            this_switch_0 = switch_to_0{tt};
            % loop through each switch and get click-triggered response
            for m=1:length(this_switch_0)
                
                % find time bin of response vector that most closely matches switch
                [t_start, t_end, ind] = get_lag_params(t,this_switch_0(m),n_lags);
                % add this STR to the matrix
                if sum(ind)>0
                    n_computed_switches_to_0 = n_computed_switches_to_0+1;
                    STR_left(n_computed_switches_to_0,max(1,t_start):t_end,:) = s2(ind,:);
                end
            end
        end
        % repeat for switches to state 1
        if ~isempty(switch_to_1{tt})
            this_switch_1 = switch_to_1{tt};
            for m=1:length(this_switch_1)
                n_computed_switches_to_1 = n_computed_switches_to_1+1;
                [t_start, t_end, ind] = get_lag_params(t,this_switch_1(m),n_lags);
                % add this STR to the matrix
                if sum(ind)>0
                    STR_right(n_computed_switches_to_1,max(1,t_start):t_end,:) = s2(ind,:);
                end
                
                
            end
        end
    end
end

% populate output structure
dprime              = (nanmean(STR_right)-nanmean(STR_left)) ./ ...
    (.5 .* (nanvar(STR_right)+nanvar(STR_left)));

res.STR_right_real  = STR_right(:,:,1);
res.STR_left_real   = STR_left(:,:,1);
res.STR_right_shuff = STR_right(:,:,2:end);
res.STR_left_shuff  = STR_left(:,:,2:end);
res.dprime_real     = squeeze(dprime(1,:,1))';
res.dprime_shuff    = squeeze(dprime(1,:,2:end));
res.good_lags       = ~isnan(sum(res.dprime_shuff,2));
res.pval            = nan(size(res.dprime_real));
res.pval(res.good_lags) = nanmean(res.dprime_real(res.good_lags)'  > res.dprime_shuff(res.good_lags,:)');
%res.pval(res.good_lags) = nanmean(abs(res.dprime_real(res.good_lags))'  > abs(res.dprime_shuff(res.good_lags,:))');

res.lags            = lags;
res.params          = p;
res.cellid          = data.cellid;

% save results
if ~isempty(p.save_dir) && p.save_file
    save(save_path,'res');
end



function [t_start t_end ind] = get_lag_params(t,this_switch_t,n_lags)
dt = t(2)-t(1);
n_timebins = length(t);
this_t = t(1):dt:(this_switch_t+dt);

%find time bin of response vector that most closely matches switch
[~, switch_time_ind] = min(abs(this_t - this_switch_t + dt*0.5));
t_start = n_timebins-switch_time_ind+1;
t_end   = n_lags - switch_time_ind + 1;
ind = (t_start:t_end)>0;






