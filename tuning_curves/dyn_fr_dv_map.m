function res = dyn_fr_dv_map(cellid,varargin)
% function res = dyn_fr_dv_map(cellid, varargin)
p = inputParser;
addParameter(p,'lag',        0);
addParameter(p,'frbins',     100);
addParameter(p,'use_nans',   0);
addParameter(p,'mask_after_stim',   true);
addParameter(p,'mask_during_stim',  false);
addParameter(p,'average_a_vals',    true)
addParameter(p,'average_fr',        true);
addParameter(p,'end_mask_s',        0);
addParameter(p,'dt',        0.05);
addParameter(p,'max_t',     1.9);
addParameter(p,'t0s',    []);
addParameter(p,'alignment', 'stimstart');   % string matching one of align_strs in dyn_align_LUT
addParameter(p,'trialnums',  []);
addParameter(p,'krn_width',  []);      % forces neural data to have this kernel width; if empty, uses whatever is in data    'krn_type'   []          );      % forces neural data to have this kernel type; if empty, uses whatever is in data
addParameter(p,'krn_type',   'halfgauss');      % forces neural data to have this kernel type
addParameter(p,'fr_dt',      0.0250);      % forces neural data to have this bin size; if empty, uses whatever is in data
addParameter(p,'norm_type',     'none');      % type of firing rate normalization; 'div' is divisive, 'z' is z-score, 'none' no normalization
addParameter(p,'frates',   []);
addParameter(p,'ft',   []);
addParameter(p,'shuffle_trials',   0);
addParameter(p,'shuffle_trials_fixing_choice',   0);
addParameter(p,'shuffle_switches',   0);
addParameter(p,'model',   []);
addParameter(p,'n_dv_bins',   []);
addParameter(p,'which_switch',   []);
addParameter(p,'which_trials',   []);
addParameter(p,'data',   []);
addParameter(p,'do_svd',   1);
addParameter(p,'demean_frates',   0);
addParameter(p,'zscore_frates',   0);
addParameter(p,'hits_only',   0);
addParameter(p,'errs_only', 0);
addParameter(p,'min_pre_dur',   0);
addParameter(p,'min_post_dur', 0);
addParameter(p,'min_switch_t',   0);
addParameter(p,'max_switch_t', 2);
addParameter(p,'verbose', 0);



parse(p,varargin{:});
p = p.Results;

if p.hits_only+p.errs_only == 2
    error('can''t do hits and errs only')
end

if length(cellid) > 1
    for cc = 1:length(cellid)
        res(cc) = dyn_fr_dv_map(cellid(cc),varargin{:});
    end
    return
end

if p.mask_during_stim && p.mask_after_stim
    error('You cant mask stimulus firing and after-stimulus firing')
end

dp = set_dyn_path;

% find alignment index for firing rates
[align_strs, align_args] = dyn_align_LUT;
align_ind = strmatch(p.alignment,align_strs,'exact');

% set up time axis for binning firing rates
if isempty(p.t0s)
    p.t0s     = p.dt:p.dt:p.max_t - p.lag;
end
if length(p.frbins) == length(p.t0s)
    error('p.t0s are same length as frbins, could get confusing')
end

% pull the data for this cell recording + behavior session
if isempty(p.data)
    data    = dyn_cell_packager(cellid);
else
    data    = p.data;
end

if ~isfield(data.trials, 'trialnums')
    if all(isempty([p.krn_width, p.krn_type, p.fr_dt]))
        data = dyn_cell_packager(cellid,'repack',1);
    else
        data = dyn_cell_packager(cellid, 'krn_width', p.krn_width, ...
            'krn_type', p.krn_type, 'bin_size', p.fr_dt,'repack',1);
    end
end

% pull time axis associated with firing rates
if isempty(p.ft)
    ft      = data.frate_t{align_ind};
else 
    ft = p.ft;
end

if isempty(p.trialnums)
    trialnums = data.trials.trialnums;
else
    trialnums = p.trialnums;
end
if isempty(p.which_trials)
    which_trials = true(size(data.trials.trialnums));
else
    which_trials = p.which_trials;
end
if p.hits_only
    which_trials = which_trials & data.trials.hit(which_trials);
end
if p.errs_only
    which_trials = which_trials & ~data.trials.hit(which_trials);
end
trialnums = trialnums(which_trials);



% pull firing rates from the appropriate alignment
if isempty(p.frates)
    frates = data.frate{align_ind}(which_trials,:);
else
    frates = p.frates(which_trials,:);
end

T       = data.trials.T(which_trials);

% normalize firing rates
if isempty(p.norm_type)
    p.norm_type = 'none';
end

switch p.norm_type
    case 'none'  
        fr_norm = frates;
    case 'div'
        fr_norm = frates ./ data.norm_mean;
    case 'z'
        fr_norm = (frates - data.norm_mean) ./ data.norm_std;
    otherwise
        error('Do not recognize norm_type')
end

if p.demean_frates | p.zscore_frates
    fr_norm  = fr_norm - nanmean(fr_norm);
end
if p.zscore_frates
    fr_norm  = fr_norm ./ nanstd(frates);
end

if isscalar(p.frbins)
    nfrbins = p.frbins;
    p.frbins = [];
elseif isempty(p.frbins)
    nfrbins = 20;
end

if isempty(p.frbins)
    p.frbins = linspace(min(fr_norm(:)), max(fr_norm(:)), nfrbins);
end

% find time offset between model time and fr time
% (model is stim_start aligned, fr alignment is chosen by user)
ref_ev_ind  = find(cellfun(@(x) strcmp('ref_event',x), align_args{align_ind})) + 1;
ref_event   = align_args{align_ind}{ref_ev_ind};

ref_times   = data.trials.(ref_event)(which_trials);
ss_times    = data.trials.stim_start(which_trials); % stim_start is ref event for model
mt_offset   = ref_times - ss_times;  % offset between model and fr_t on each trial

sessid      = data.sessid;

if ~isempty(p.which_switch)
    clear_bad_strengths = 1;
    bad_strength    = 0;
    fit_line = 1;
    exclude_final = 0;
    final_only = 0;
    [switch_to_0, switch_to_1, ~, vd] = ...
        get_switches(cellid, 'which_switch',p.which_switch,...
        'clear_bad_strengths', clear_bad_strengths, ...
        'bad_strength', bad_strength, 'fit_line', fit_line,...
        'exclude_final', exclude_final, 'final_only', final_only,...
        'min_pre_dur',p.min_pre_dur,'min_post_dur',p.min_post_dur,...
        'which_trials',which_trials,'min_switch_t',p.min_switch_t,...
        'max_switch_t',p.max_switch_t);
    
        
else
    switch_to_0 = [];
    switch_to_1 = [];
end
assert(length(data.trials.gamma) == length(switch_to_0));


if isempty(p.model)
    % only grab the trials that are in the data struct
    model = get_data_model_p(data, trialnums);
else
    %warning('temporary fix to choose model trials')
    if length(p.model) > length(which_trials)
        model = get_data_model_p(p.model, trialnums(which_trials));
    elseif length(p.model) == length(which_trials)
        model = get_data_model_p(p.model, which_trials);
    else
        error('not sure what to do with model trialnums')
    end
    p.model = [];
end
    
modelT  = cellfun(@(x) x.T(end), {model.posterior});

if ~(all(abs(modelT(:)-T(:)) < .1))
    warning('model T and T don''t match up perfectly')
end
%{
there's some discrepancy between the time in the model posterior and the
trial duration T which corresponds. you can see that here
spike_data_dir = dp.spikes_dir;
[~, vec_data, ~, ~] = get_behavior_data(spike_data_dir, cellid, sessid);
S       = load_session_data(sessid);
Sdata   = format_data(S);
%}

if p.shuffle_trials | p.shuffle_trials_fixing_choice | p.shuffle_switches
    Ttiles = [min(T); percentile(T,25:25:75); max(T)+eps];
    if p.shuffle_trials_fixing_choice
        choice_group = data.trials.rat_dir;
    else
        choice_group = ones(size(data.trials.rat_dir));
    end
    unique_choice_group = unique(choice_group);
    for gg = 1:length(unique_choice_group)
        this_group = choice_group == unique_choice_group(gg);
        % resort wp.ithin duration quintiles
        for tt = 1:length(Ttiles)-1
            goodT   =  T >= Ttiles(tt) & T < Ttiles(tt+1);
            this_ind = find(this_group & goodT);
            [~, neworder]       = sort(rand(length(this_ind),1));
            if p.shuffle_trials
                fr_norm(this_ind,:) = fr_norm(this_ind(neworder),:);
            end
            if p.shuffle_switches
                switch_to_0(this_ind) = switch_to_0(neworder);
                switch_to_1(this_ind) = switch_to_1(neworder);
            end
        end
    end
end


% compute the joint probability distribution p(t,f,a)
res = compute_joint_frdvt(fr_norm, ft, model, p.t0s, mt_offset, ...
    'frbins', p.frbins, 'use_nans', p.use_nans, 'lag', p.lag,...
    'mask_after_stim', p.mask_after_stim, 'mask_during_stim', p.mask_during_stim,...
    'average_a_vals', p.average_a_vals, 'average_fr', p.average_fr,...
    'end_mask_s', p.end_mask_s, ...
    'switch_to_0', switch_to_0, 'switch_to_1',switch_to_1,...
    'verbose', p.verbose);
%
if ~isempty(p.n_dv_bins)
    res = rebin_joint(res, p.n_dv_bins, 'use_nans', p.use_nans);
    
end
% compute the 2d tuning, i.e. E[f(t,a)] = p(t,a,f)*f
res = compute_tuning_from_joint(res);
% compute average 1d tuning, i.e. E[f(a)] = <E[f(t,a)]>
res = compute_dv_tuning(res);
% run svd and et rank 1 matrix approximation of E[f(t,a)]
if p.do_svd
    res = compute_rank1_fgta_approx(res);
end

res.cellid = cellid;
res.which_switch = p.which_switch;
res.params = p;


function res = compute_joint_frdvt(fr_norm, ft, model, t0s, mt_offset, varargin)
p = inputParser;
addParameter(p, 'frbins', 0:.5:100);
addParameter(p, 'use_nans', 0);
addParameter(p, 'lag', 0);
addParameter(p, 'mask_after_stim', true);
addParameter(p, 'mask_during_stim', false);
addParameter(p, 'average_a_vals', true)
addParameter(p, 'average_fr', true);
addParameter(p, 'end_mask_s', 0);
addParameter(p, 'switch_to_0', '');
addParameter(p, 'switch_to_1', '');
addParameter(p, 'verbose', 0);
parse(p,varargin{:});
p = p.Results;
frbins = p.frbins;

end_mask_s  = p.end_mask_s;
switch_to_0 = p.switch_to_0;
switch_to_1 = p.switch_to_1;

use_switches = ~isempty(switch_to_0) & ~isempty(switch_to_1);

if use_switches
    n_switch_to_0 = cellfun(@length, switch_to_0);
    n_switch_to_1 = cellfun(@length, switch_to_1);
    has_switch  = n_switch_to_0 + n_switch_to_1;
    includes    = has_switch(:);
end

assert(size(fr_norm,1)==length(model))
assert(size(fr_norm,1)==length(mt_offset) | length(mt_offset) == 1)

% check that the model posterior is binned the same for each trial
[model, constant_a, allsame] = get_data_model_p(model);
if ~allsame
    numa = length(constant_a);
else
    numa = numel(model(1).posterior.avals);
end

abin_width  = diff(constant_a(1:2));
mass_t      = zeros(numel(t0s),1);
mass_tfa    = zeros(numel(t0s), numel(frbins), numa);
if p.use_nans
    mass_tfa    = mass_tfa*nan;
end;

t0_durs     = [diff(t0s) 0];
last_durs   = [0 t0_durs];
counter     = 0;
if p.verbose, fprintf('starting trial 1...'); end
for ti = 1:numel(model)
    if mod(ti,100)==0 & p.verbose
        fprintf('%i...',ti)
    end
    % get this trials model predictions and firing rates
    a   = model(ti).posterior.avals;
    t   = model(ti).posterior.T - mt_offset(ti);
    pr  = model(ti).posterior.pdf;
    fr  = fr_norm(ti,:);
    
    if use_switches
        this_switches   = sort([switch_to_0{ti} switch_to_1{ti}]);
        this_nswitches  = length(this_switches);
    else
        this_switches   = 0;
        this_nswitches  = 1;
    end
    % get last time point where we have a firing rate for this trial
    last_ft = ft(find(~isnan(fr),1,'last')) - end_mask_s;
    
    for si = 1:this_nswitches
        counter     = counter + 1;
        switch_t0s  = this_switches(si) + t0s;
        
        pre_mask    = ft(1);
        post_mask   = last_ft;
        
        if si < this_nswitches
            post_mask  = this_switches(si+1);
        end
        
        if si > 1
            pre_mask  = this_switches(si-1);
        end
        
        good_tinds = find(switch_t0s >= pre_mask & switch_t0s < post_mask);
        for tx = good_tinds
            t0          = switch_t0s(tx);
            last_dur    = 0;
            this_dur    = 0;
            if tx > good_tinds(1)
                last_dur = last_durs(tx);
            end
            if tx < good_tinds(end)
                this_dur = t0_durs(tx);
            end
            if t0 + this_dur > last_ft
                continue
            end
            
            % figure out which indices to include for posterior and firing rate
            if (t0 + this_dur/2) > last_ft
                this_dur_nolag = 2*(last_ft - t0);
            else
                this_dur_nolag = this_dur;
            end
            if (t0 + this_dur/2 + p.lag) > last_ft
                this_dur_lag = 2*(last_ft - t0 - p.lag);
            else
                this_dur_lag = this_dur;
            end
            
            t_start     = (t0 - last_dur/2);
            model_t_end = (t0 + this_dur_nolag/2);
            fr_t_end    = (t0 + this_dur_lag/2);
            % Determine relevant timepoints for accumulator value a
            if p.average_a_vals
                model_t_idx = t > t_start  & t < model_t_end;
            else
                [~, model_t_idx]   = min(abs(t-t0));
            end
            
            % Get firing rate r0 for this time point
            if p.average_fr
                %warning(t0~=fr_t_end);
                % average of interpolated values at 3 timepoints
                r0 = nanmean(interp1(ft, fr, p.lag + [t_start, t0, fr_t_end]));% The mean firing rate during the time window
            else
                % interpolated value at timepoint
                r0 = interp1(ft, fr, t0 + p.lag);
            end
            % get index into joint distribution for this firing rate
            [~, fbin] = min(abs(frbins-r0));
            
            % decide whether to include this timepoint
            if p.mask_after_stim
                timepoint = t(end)-(t0+p.lag);%-mt_offset(ti);
                include_fr = ~isnan(r0) & (timepoint > 0);
            elseif p.mask_during_stim
                timepoint = t(end)-(t0+p.lag);%)-mt_offset(ti);
                include_fr = ~isnan(r0) & ( timepoint< 0);
            else
                include_fr = ~isnan(r0);
            end
            
            % decide whether to add mass to joint distribution for this time bin and firing rate
            if include_fr
                % select middle numx timepoints if p is different size
                thisx   = size(pr,2);
                sd      = ceil(thisx/2) - floor(numa/2);
                ed      = ceil(thisx/2) + floor(numa/2);
                this_pr = mean(pr(model_t_idx, sd:ed),1);
                % normalize p(a) for this time point so integrates to 1
                this_pr = this_pr./sum(this_pr)/abin_width;
                % add p(a) mass to this time and firing rate bin
                if p.use_nans && isnan(mass_tfa(tx, fbin, 2))  % this should only happen is use_nans was 1
                    mass_tfa(tx, fbin, :) = this_pr;
                else
                    mass_tfa(tx, fbin, :) = squeeze(mass_tfa(tx, fbin,:))' ...
                        + this_pr;
                end
                % keep track of how many trials contribute to each time bin
                mass_t(tx)  = mass_t(tx) + 1;
            end
        end
    end
end
res.numa        = numa;
res.t0s         = t0s;
res.frbins      = p.frbins;
res.dv_axis     = constant_a;
res.mass_tfa    = mass_tfa;
res.mass_t      = mass_t;
res.pjoint_tfa  = mass_tfa ./ mass_t;




function res = compute_tuning_from_joint(res)
%% get joint and conditional distributions
pj_given_a      = zeros(numel(res.t0s), numel(res.frbins), res.numa);
pj_given_fr     = zeros(numel(res.t0s), numel(res.frbins), res.numa);
fr_given_ta     = zeros(numel(res.t0s), res.numa);
fr_var_given_ta = zeros(numel(res.t0s), res.numa);
a_given_tfr     = zeros(numel(res.t0s), numel(res.frbins));

for tx = 1:numel(res.t0s) % iterate over timepoints
    
    this_pfa = squeeze(res.pjoint_tfa(tx,:,:));
    % normalize Pjoint to get firing rate tuning curve wrt a
    this_fnorm = ones(size(this_pfa,1),1)*sum(this_pfa,1);
    pj_given_a(tx,:,:) = this_pfa ./ this_fnorm;
    % normalize Pjoint to get a distribution wrt firing rate
    this_anorm = sum(this_pfa,2)*ones(1,size(this_pfa,2));
    pj_given_fr(tx,:,:)     = this_pfa ./ this_anorm;
    
    fr_given_ta(tx,:)        = squeeze(pj_given_a(tx,:,:))'*res.frbins';
    fr_var_given_ta(tx,:)   = (squeeze(pj_given_a(tx,:,:))'*(res.frbins'.^2)) - ...
        (fr_given_ta(tx,:)'.^2);
    a_given_tfr(tx,:)       = squeeze(pj_given_fr(tx,:,:))*res.dv_axis';
end;

res.pj_given_a      = pj_given_a;
res.pj_given_fr     = pj_given_fr;
res.fr_given_ta      = fr_given_ta;
res.fr_var_given_ta  = fr_var_given_ta;
res.a_given_tfr      = a_given_tfr;


function res = compute_dv_tuning(res,varargin)
p = inputParser;
addParameter(p, 'min_fr_mod', 0)
parse(p,varargin{:});
ops = p.Results;

ntimes      = size(res.fr_given_ta,1);
time_bins   = 1:ntimes;
fgta        = res.fr_given_ta;

% compute 1d tuning curve by subtracting time average from fgta
% normalizing, and then averaging over time

% subtract mean tuning curve from each timepoint
fgta_resid  = fgta - nanmean(fgta,2);
fgta_resid_n   = fgta_resid - min(fgta_resid,[],2);
fgta_rn_max    = max(fgta_resid_n, [], 2);
fgta_resid_n   = fgta_resid_n ./ fgta_rn_max;
% compute range of firing rates at each timepoint
frm_time    = max(fgta_resid,[],2)-min(fgta_resid,[],2);
% decide which timepoints to include
good_tind   = frm_time > ops.min_fr_mod;
% average residual normalized fgta over time
fga_rn_tmn  = nanmean(fgta_resid_n(good_tind,:),1);
% average residual normalized fgta over time
fga_rn_std = nanstderr(fgta_resid_n(good_tind,:),1);
fga_rn_std = fga_rn_std ./ fgta_rn_max;

% compute 1d tuning curve by averaging fgta over time
fga_tmn      = nanmean(fgta(time_bins,:) - ...
    nanmean(fgta(time_bins,:),2),1);
fga_std  = nanstderr(fgta(time_bins,:) - ...
    nanmean(fgta(time_bins,:),2),1);
fr_mod   = max(fga_tmn) - min(fga_tmn);
fga_tmn_n       = fga_tmn -  min(fga_tmn);
fga_tmn_n_max   = max(fga_tmn_n);
fga_tmn_n       = fga_tmn_n / fga_tmn_n_max;
fga_std_n       = fga_std ./ fga_tmn_n_max;

fga_resid_tmn   = nanmean(fgta_resid(time_bins,:));
fga_resid_std  = nanstderr(fgta_resid(time_bins,:));
fga_tmn         = nanmean(fgta(time_bins,:));
fga_std        = nanstderr(fgta(time_bins,:));

res.fgta_resid      = fgta_resid;
res.fgta_resid_n    = fgta_resid_n;
res.fga_rn_tmn      = fga_rn_tmn;
res.fga_rn_std      = fga_rn_std;
res.fga_resid_tmn   = fga_resid_tmn;
res.fga_resid_std   = fga_resid_std;
res.fga_tmn         = fga_tmn;
res.fga_std         = fga_std;
res.fga_tmn_n       = fga_tmn_n;
res.fga_std_n       = fga_std_n;
res.good_tind       = good_tind;
res.frm_time        = frm_time;
res.fr_mod          = fr_mod;




function res = rebin_joint(res, n_dv_bins,varargin)
p = inputParser;
addParameter(p, 'use_nans', 0);
parse(p,varargin{:});
p = p.Results;


x   = res.dv_axis;
t0s = res.t0s;
frbins = res.frbins;
pj_fine     = res.pjoint_tfa;
ma_fine    = res.mass_tfa;

if n_dv_bins > size(pj_fine,3)
    error('You are binning accumulation values more finely than the model generates')
end

% set up new dv bins

if length(n_dv_bins) > 1
    dv_bin_edges = n_dv_bins;
    dvbe = dv_bin_edges;
    dvax = [ dvbe(1:end-1) + 0.5*(diff(dvbe)) ];
    dv_axis = [min(x)+(dvax(1)-min(x))/2 dvax max(x)+(dvax(end)-max(x))/2];
else
x_bound_margin = (x(2)-x(1)) * 0.5;
    dv_bin_edges = linspace(min(x)+x_bound_margin,max(x)-x_bound_margin,n_dv_bins-1);

    dv_axis = dv_bin_edges(1:end-1) + 0.5*(diff(dv_bin_edges));
    dv_axis = [min(x) dv_axis max(x)];  % add bound bins for axis
end

n_dv_bins = length(dv_axis);

new_pj = zeros(size(pj_fine,1), size(pj_fine,2), length(dv_axis));
new_ma = zeros(size(pj_fine,1), size(pj_fine,2), length(dv_axis));
if p.use_nans
    new_pj = nan * new_pj;
    new_ma = nan * new_pa;
end
pj_given_ta     = zeros(size(pj_fine,1), size(pj_fine,2), length(dv_axis));
fr_given_ta     = zeros(numel(t0s), numel(dv_axis));
fr_var_given_ta = zeros(numel(t0s), numel(dv_axis));

for time_i=1:numel(t0s)   
    % Deal with the first bin
    first_dv_bins       = x < dv_bin_edges(1);
    if sum(first_dv_bins)==1
        new_pj(time_i,:,1) = pj_fine(time_i,:,first_dv_bins);
        new_pj(time_i,:,1) = pj_fine(time_i,:,first_dv_bins);
        new_ma(time_i,:,1) = ma_fine(time_i,:,first_dv_bins);
        new_ma(time_i,:,1) = ma_fine(time_i,:,first_dv_bins);
    else
        new_pj(time_i,:,1) = nansum(pj_fine(time_i,:,first_dv_bins),3);
        new_ma(time_i,:,1) = nansum(ma_fine(time_i,:,first_dv_bins),3);
    end
    edge_match          = x == dv_bin_edges(1);
    new_pj(time_i,:,1) = new_pj(time_i,:,1) + 0.5*nansum(pj_fine(time_i,:,edge_match),3);
    new_ma(time_i,:,1) = new_ma(time_i,:,1) + 0.5*nansum(ma_fine(time_i,:,edge_match),3);

    % Deal with middle bins
    for j=2:n_dv_bins-1
        these_dv_bins       = x>dv_bin_edges(j-1) & x<dv_bin_edges(j);
        if sum(these_dv_bins) == 1
            new_pj(time_i,:,j) = pj_fine(time_i,:,these_dv_bins); % Pjoints is t0 x FR x dv
            new_ma(time_i,:,j) = ma_fine(time_i,:,these_dv_bins); % Pjoints is t0 x FR x dv

        else
            new_pj(time_i,:,j) = nansum(pj_fine(time_i,:,these_dv_bins),3); % Pjoints is t0 x FR x dv
            new_ma(time_i,:,j) = nansum(ma_fine(time_i,:,these_dv_bins),3); % Pjoints is t0 x FR x dv
        end
        % deal with bin edges that exactly match DV centers; for
        % these case, split joint dist in half between adjacent bins
        edge_match          = x==dv_bin_edges(j-1) | x==dv_bin_edges(j);
        new_pj(time_i,:,j) = new_pj(time_i,:,j) + 0.5*nansum(pj_fine(time_i,:,edge_match),3);
        new_ma(time_i,:,j) = new_ma(time_i,:,j) + 0.5*nansum(ma_fine(time_i,:,edge_match),3);
    end
    % Deal with last bin
    final_dv_bins           = x > dv_bin_edges(end);
    if sum(final_dv_bins) == 1
        new_pj(time_i,:,end)   = pj_fine(time_i,:,final_dv_bins);
        new_ma(time_i,:,end)   = ma_fine(time_i,:,final_dv_bins);
    else
        new_pj(time_i,:,end)   = nansum(pj_fine(time_i,:,final_dv_bins),3);
        new_ma(time_i,:,end)   = nansum(ma_fine(time_i,:,final_dv_bins),3);
    end
    edge_match              = x == dv_bin_edges(end);
    new_pj(time_i,:,end)   = new_pj(time_i,:,end) + 0.5*nansum(pj_fine(time_i,:,edge_match),3);
    new_ma(time_i,:,end)   = new_ma(time_i,:,end) + 0.5*nansum(ma_fine(time_i,:,edge_match),3);

    p_mass_diff = sum(pj_fine(time_i,:,:),3) - sum(new_pj(time_i,:,:),3);
    diff_mag    = sum(abs(p_mass_diff));
    if diff_mag > 1e-6
        keyboard
        figure(111); clf
  
        s(1)= subplot(211)
        temp = squeeze(new_pj(time_i,1,:));
        plot(dv_axis,temp,'.-')
        s(2)= subplot(212)
        temp = squeeze(pj_fine(time_i,1,:));
        plot(x,temp,'.-')
        linkaxes(s)
        %%
        warning('Probability mass is leaking during binning')
    end
    myPj = squeeze(new_pj(time_i,:,:));
    this_pj_given_ta = myPj ./ (ones(size(myPj,1),1)*sum(myPj,1));
    fr_given_ta(time_i,:) = this_pj_given_ta'*frbins';
    fr_var_given_ta(time_i,:) = this_pj_given_ta'*(frbins'.^2) - (fr_given_ta(time_i,:)'.^2);
    pj_given_ta(time_i,:,:)   = this_pj_given_ta;
end;


% save binned dv_axis as x
res.mass_tfa        = new_ma;
res.pjoint_tfa      = new_pj;
res.pj_given_ta     = pj_given_ta;
res.fr_given_ta     = fr_given_ta;
res.fr_var_given_ta = fr_var_given_ta;
res.dv_axis         = dv_axis;
res.numa            = numel(dv_axis);

%%
function results = dyn_fr_dv_map_(cellids, t0s, n_dv_bins, ops, varargin)
% function fr_dv_map(cellids, t0s, n_dv_bins, varargin)
%
% Function that performs analysis of firing rate's relationship to
% behavioral decision variable. It generates two figures: 1) Firing rate versus time
% parameterized on DV value in color. 2) The FR-DV relationship averaged
% across specified times
%
% Input:
%       cellids:    vector of cellids to analyze
%       t0s:        vector of times to examine. These will be times into the
%                       model's evolution relative to alignment.
%       n_dv_bins:  number of bins to use for decision variable. Bounds are
%                   included as outside bins, so this must be >=3.
%
%
% Output:
%       NONE
%
%
%
% %%% Default Parameters (change with varargin)

p = inputParser;
addParameter(p, 'lag', .2) % click lag: 0.2
addParameter(p, 'time_edges', []) % first and last time over which to calculate average; if empty use all t0s
addParameter(p, 'var_weight', true); % weight mean by inverse variance if true
addParameter(p, 'frbins', 0:0.1:50);  % bins for FR in joint distribution
addParameter(p, 'dt',  0.01); % model dt
addParameter(p, 'alignment',  'stimstart-cout-mask'); % string matching one of align_strs in dyn_align_LUT
addParameter(p, 'direction',  'backward'); % 'forward' or 'backward'
addParameter(p, 'krn_width',  []); % forces neural data to have this kernel width; if empty, uses whatever is in data
addParameter(p, 'fr_dt',  0.0250); % forces neural data to have this bin size; if empty, uses whatever is in data
addParameter(p, 'krn_type',  'halfgauss'); % forces neural data to have this kernel type
addParameter(p, 'norm_type',  'div'); % type of firing rate normalization; 'div' is divisive, 'z' is z-score
addParameter(p, 'flip_bool',  '~data.prefsideright{align_ind}'); % boolean choice for each cellid of whether to flip sign of DV (string evaluated in workspace); default makes pref direction the higher DV value
addParameter(p, 'force_frdv',  0); % force extraction from database rather than using saved file
addParameter(p, 'force_bin',  0);
addParameter(p, 'force_dv',  0);
addParameter(p, 'save_path',  '');
addParameter(p, 'save_filename',  'unknown'); % a name to specify this particular map (perhaps the group of cells used)
addParameter(p, 'save_map',  true); % whether to save the map as one big file
addParameter(p, 'axes_in',  []); % axis to use for time plot; if empty, create new figure
addParameter(p, 'axes_in_ta',  []); % axis to use for time-average plot; if empty, create new figure
addParameter(p, 'ta_color',  'k'); % time-average axis color
addParameter(p, 'ta_offset',  0); % offset in x-axis data points for time-average plot
addParameter(p, 'dv_bin_edges',  []); % if empty, makes n_dv_bins quantilized bins; otherwise, override with this
addParameter(p, 'norm_ta',  true); % whether to "normalize" time average between zero and one
addParameter(p, 'ta_marker',  ''); % marker type for time average plot
addParameter(p, 'connect_bound',  true); % whether to connect the bound bin on plot to non-bound DVs with dotted line
addParameter(p, 'plot_errorbars',  true); % plot error bars on time plot
addParameter(p, 'ta_use_colors',  true); % whether to use different colors for different ta points
addParameter(p, 'plot_limit',  100); % most extreme DV value to plot (usually make this correspond to lowest bound for any rat)
addParameter(p, 'Markersize',  2); % markersize of plot
addParameter(p, 'bound',  inf); % bound for slope-fitting purposes
addParameter(p, 'n_iter',  1); % number of iterations for refining estimate of DV
addParameter(p, 'param_scale_num',  1); % parameter number to scale
addParameter(p, 'param_scale_factor',  1); % multiplicative factor of that parameter
addParameter(p, 'do_flip', 0);
addParameter(p, 'shuffle_trials', 0);
addParameter(p, 'frates', []);
parse(p, varargin{:});
%bootstrap = false;                      % if true, sample with replacement across trials; used for bootstrap stats
struct2vars(p.Results);

if isempty(save_path)  %#ok<NODEF>
    dp          = set_dyn_path;
    save_path   = dp.celldat_dir;
    datadir     = fullfile(save_path);
end

if n_dv_bins<3
    error('You must have at least three DV bins. 1 for each bound and 1 for non-bound values.')
end

% convert time edges to bin numbers
if isempty(time_edges)
    time_bins = 1:numel(t0s);
elseif length(time_edges)~=2
    error('time_edges must be empty of a vector of length 2.')
else
    [~,first_bin]   = min(abs(t0s-time_edges(1)));
    [~,last_bin]    = min(abs(t0s-time_edges(2)));
    time_bins       = first_bin:last_bin;
end

if ~isempty(dv_bin_edges)
    n_dv_bins = length(dv_bin_edges) + 1;
end


if ~force_frdv && exist([save_path filesep save_filename '_frdvmapdata.mat'], 'file')
    load([save_path filesep save_filename '_frdvmapdata.mat']);
else
    % setup space
    n_cells             = numel(cellids);
    fga_cell            = nan(numel(t0s),n_dv_bins,n_cells);
    fvga_cell           = nan(numel(t0s),n_dv_bins,n_cells);
    fga_cell_residual   = nan(numel(t0s),n_dv_bins,n_cells);
    fga_cell_nresidual  = nan(numel(t0s),n_dv_bins,n_cells);
    fga_ta_cell         = nan(n_dv_bins,n_cells);
    fga_nta_cell        = nan(n_dv_bins,n_cells);
    fga_ta_unorm_cell   = nan(n_dv_bins,n_cells);
    fga_std_ta_cell     = nan(n_dv_bins,n_cells);
    fga_std_nta_cell    = nan(n_dv_bins,n_cells);
    fga_aa_cell         = nan(numel(t0s),n_cells);
    x_cell              = nan(n_dv_bins,n_cells);
    fr_modulation       = nan(n_cells,1);
    frm_time            = nan(numel(t0s), n_cells);
    flipdex             = zeros(n_cells,1);
    rank1_fga_resid         = nan(numel(t0s), n_dv_bins,n_cells);
    rank1_mt_n           = nan(numel(t0s), n_cells);
    rank1_ra_n           = nan(n_dv_bins,n_cells);
    rank1_variance      = nan(n_cells,1);
    
    % loop through cells
    for ci=1:n_cells
        fprintf('Joint FR-DV analysis on cell %d of %d \n', ci, numel(cellids));
        
        %        rat = bdata('select ratname from cells where cellid="{S}"',cellids(ci));
        
        [x, frbins, Pjoints, fr_given_as, fr_var_given_as] = ...
            dyn_compile_binned_database(cellids(ci), t0s, n_dv_bins, ops,...
            'lag', lag, 'krn_width', krn_width, 'fr_dt', fr_dt, 'dt', dt, 'alignment', alignment,...
            'direction', direction, 'frbins', frbins, 'krn_type', krn_type, 'norm_type', norm_type,...
            'n_iter',n_iter, 'param_scale_num', param_scale_num, ...
            'param_scale_factor', param_scale_factor, 'force_bin', force_bin,'force_dv',force_dv,...
            'datadir',datadir,'shuffle_trials',shuffle_trials,'frates',frates);%, 'bootstrap', bootstrap);
        try
            
            % Do we need to flip cell?
            fga_ta_temp = nanmean(fr_given_as(time_bins,:)- nanmean(fr_given_as(time_bins,:),2),1);
            flip_cond  = mean(fga_ta_temp(x>0)) < mean(fga_ta_temp(x<0));
            if do_flip & flip_cond
                fr_given_as     = flipdim(fr_given_as,2);
                fr_var_given_as = flipdim(fr_var_given_as,2);
                flipdex(ci)     = 1;
            end
            
            % Compute fga, and fga_residual
            fga_cell(:,:,ci)    = fr_given_as;
            fga_cell_residual(:,:,ci) = fr_given_as - nanmean(fr_given_as,2);
            fga_aa_cell(:,ci)   = nanmean(fr_given_as,2);
            fvga_cell(:,:,ci)   = fr_var_given_as;
            x_cell(:,ci)        = x;
            
            % normalize fga_residual for each time bin
            for t = 1:size(fga_cell_residual,1);
                fga_cell_nresidual(t,:,ci) = fga_cell_residual(t,:,ci);
                fga_cell_nresidual(t,:,ci) = fga_cell_nresidual(t,:,ci) - min(fga_cell_nresidual(t,:,ci));
                fga_cell_nresidual(t,:,ci) = fga_cell_nresidual(t,:,ci)./max(fga_cell_nresidual(t,:,ci));
            end
            frm_time(:,ci) = max(fga_cell_residual(:,:,ci)') - min(fga_cell_residual(:,:,ci)');
            fga_nta     = nanmean(fga_cell_nresidual(time_bins(frm_time(:,ci) > ops.min_fr_mod),:,ci),1);
            fga_std_nta = nanstderr(fga_cell_nresidual(time_bins(frm_time(:,ci) > ops.min_fr_mod),:,ci),1);
            min_val = min(fga_nta);
            fga_nta = fga_nta - min_val;
            max_val = max(fga_nta);
            fga_nta = fga_nta./max_val;
            fga_std_nta = fga_std_nta ./max_val;
            
            % Find 1D tuning curve for each cell
            fga_ta      = nanmean(fr_given_as(time_bins,:) - ...
                nanmean(fr_given_as(time_bins,:),2),1);
            fga_std_ta  = nanstderr(fr_given_as(time_bins,:) - ...
                nanmean(fr_given_as(time_bins,:),2),1);
            fga_ta_unorm = fga_ta;
            min_val     = min(fga_ta);
            fga_ta      = fga_ta - min_val;
            max_val     = max(fga_ta);
            fga_ta      = fga_ta ./ max_val;
            fga_std_ta  = fga_std_ta ./ max_val;
            fga_ta_cell(:,ci)       = fga_ta;
            fga_nta_cell(:,ci)      = fga_nta;
            fga_std_ta_cell(:,ci)   = fga_std_ta;
            fga_std_nta_cell(:,ci)   = fga_std_nta;
            fga_ta_unorm_cell(:,ci) = fga_ta_unorm;
            fr_modulation(ci) = max(fga_ta_unorm) - min(fga_ta_unorm);
            
            % fit sigmoid to 1D tuning curve - using direct averaging
            try
                warning('off','all')
                [betas,resid,jacob,sigma,mse] = nlinfit(x,fga_ta,@dyn_sig,[0, 1/3]);
                warning('on','all')
                betas_cell(ci,:)    = betas;
                delta               = sqrt(diag(sigma)) * tinv(1-0.05/2,sum(~isnan(x))-4);
                sigmas_cell(ci,:)   = delta;
                slope_cell(ci)      = betas(2)/4;
            catch
                betas_cell(ci,:)    = nan;
                sigmas_cell(ci,:)   = nan;
                slope_cell(ci)      = nan;
            end
            % fit sigmoid to 1D tuning curve - using normalized / time bin
            try
                warning('off','all')
                [betas,resid,jacob,sigma,mse] = nlinfit(x,fga_nta,@dyn_sig,[0, 1/3]);
                warning('on','all')
                nbetas_cell(ci,:)    = betas;
                ndelta               = sqrt(diag(sigma)) * tinv(1-0.05/2,sum(~isnan(x))-4);
                nsigmas_cell(ci,:)   = delta;
                nslope_cell(ci)      = betas(2)/4;
            catch
                nbetas_cell(ci,:)    = nan;
                nsigmas_cell(ci,:)   = nan;
                nslope_cell(ci)      = nan;
            end
            % fit sigmoid to 1D tuning curve for each normalized time bin
            if ops.fit_time_sigmoid
                warning('off','all')
                for t = 1:size(fga_cell_nresidual,1);
                    try
                        [betas,resid,jacob,sigma,mse] = nlinfit(x,fga_cell_nresidual(t,:,ci),@dyn_sig,[0, 1/3]);
                        nbetas_time(ci,t,:)    = betas;
                        ndelta                      = sqrt(diag(sigma)) * tinv(1-0.05/2,sum(~isnan(x))-4);
                        nsigmas_time(ci,t,:)   = delta;
                        nslope_time(ci,t)      = betas(2)/4;
                    catch
                        nbetas_time(ci,t,:)    = nan;
                        nsigmas_time(ci,t,:)   = nan;
                        nslope_time(ci,t)      = nan;
                    end
                end
                warning('on','all')
            end
            
            % SVD analysis
            
            [u,s,v]     = svd(fga_cell_residual(:,:,ci));
            s_squared   = diag(s).^2;
            s1(:,1)     = s(1);
            u1(:,ci)    = u(:,1);
            v1(:,ci)    = v(:,1);
            %s1(2:end)   = 0;
            rank1_fga_resid(:,:,ci) = u(:,1)*s(1)*v(:,1)';
            alpha       = 1/range(v(:,1));
            beta        = s1(1,1)/alpha;
            
            rank1_mt_n(:,ci) = u(:,1)*beta;
            rank1_ra_n(:,ci) = v(:,1)*alpha;
            for ii = 1:5
                rank_variance(ci,ii) = sum(s_squared(1:ii))./sum(s_squared);
            end
            
            tuning_alpha(ci) = alpha;
            tuning_beta(ci) = beta;
            % sign of
            if median(rank1_mt_n(:,ci)) < 0
                rank1_mt_n(:,ci) = -u(:,1)*beta;
                rank1_ra_n(:,ci) = -v(:,1)*alpha;
            end
            % fit sigmoid to rank 1 tuning curve
            try
                warning('off','all')
                svdta = rank1_ra_n(:,ci)';
                svdta = svdta - min(svdta);
                svdta = svdta./max(svdta);
                [betas,resid,jacob,sigma,mse] = nlinfit(x,svdta,@dyn_sig,[0, 1/3]);
                warning('on','all')
                svd_betas_cell(ci,:)    = betas;
                svd_delta               = sqrt(diag(sigma)) * tinv(1-0.05/2,sum(~isnan(x))-4);
                svd_sigmas_cell(ci,:)   = svd_delta;
                svd_slope_cell(ci)      = betas(2)/4;
            catch
                svd_betas_cell(ci,:)    = nan;
                svd_sigmas_cell(ci,:)   = nan;
                svd_slope_cell(ci)      = nan;
            end
            
        catch me
            disp('something crashed in dyn_fr_dv_map')
            keyboard
        end
    end
    
    results.fvga_cell           = fvga_cell;
    results.fga_cell            = fga_cell;
    results.fga_cell_residual   = fga_cell_residual;
    results.fga_cell_nresidual  = fga_cell_nresidual;
    results.fga_ta_cell         = fga_ta_cell;
    results.fga_nta_cell        = fga_nta_cell;
    results.fga_std_ta_cell     = fga_std_ta_cell;
    results.fga_std_nta_cell    = fga_std_nta_cell;
    results.fga_ta_unorm_cell   = fga_ta_unorm_cell;
    results.fga_aa_cell         = fga_aa_cell;
    results.betas_cell          = betas_cell;
    results.sigmas_cell         = sigmas_cell;
    results.slope_cell          = slope_cell;
    results.nbetas_cell         = nbetas_cell;
    results.nsigmas_cell        = nsigmas_cell;
    results.nslope_cell         = nslope_cell;
    if ops.fit_time_sigmoid
        results.nbetas_time         = nbetas_time;
        results.nsigmas_time        = nsigmas_time;
        results.nslope_time         = nslope_time;
    end
    results.fr_modulation       = fr_modulation;
    results.frm_time            = frm_time;
    results.flipdex             = flipdex;
    results.x_cell              = x_cell;
    results.tuning_cell         = rank1_fga_resid;
    results.tuning_fr           = rank1_mt_n;
    results.tuning_fa           = rank1_ra_n;
    results.rank1_variance      = rank1_variance;
    results.svd_betas_cell      = svd_betas_cell;
    results.svd_sigmas_cell     = svd_sigmas_cell;
    results.svd_slope_cell      = svd_slope_cell;
    results.tuning_alpha        = tuning_alpha;
    results.tuning_beta         = tuning_beta;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finished iterating across cells, now population analysis
ops.var_weight        = var_weight;
ops.cell_dex          = results.fr_modulation > ops.min_fr_mod;
results.cell_dex    = ops.cell_dex;
results.cell_num    = find(ops.cell_dex);
results.t0s         = t0s;
results.dv_axis     = mean(results.x_cell,2);
results.time_bins   = time_bins;
results.cellid      = cellids;

try
    results             = population_analysis(results,ops);
catch me
    disp(me.message);
end
% plot population summary figures, only plots good cells
try
    plot_population(results);
catch me
    disp(me.message);
end

numgood = sum(ops.cell_dex);
fracgood = 100*(numgood/numel(cellids));
disp(['Had ' num2str(numgood) ' good cells out of ' num2str(numel(cellids)) ' : ' num2str(fracgood) '%'])

% Save
if save_map,
    if ~exist(save_path, 'dir'), mkdir(save_path); end;
    save([save_path filesep save_filename '_frdvmapdata.mat'], 'results');
end



