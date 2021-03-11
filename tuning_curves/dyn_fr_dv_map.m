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
    assert(length(data.trials.gamma) == length(switch_to_0));
        
else
    switch_to_0 = [];
    switch_to_1 = [];
end



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



