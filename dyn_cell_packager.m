% [data] = cell_packager(cellid, {'force', 0}, {'save_package', 1})
%
% Loads basic data for a cellid, including spike rate functions for
% different alignments. This is one processing step beyond
% package_pbups_phys function. If a saved file exists, uses it.
%
% OPTIONAL PARAMS:
% ----------------
%
%   repack    Default 0. If 0 and saved file exists, just load save file,
%             if 1 recompile the data.
%
%   save_package   Default 1.  Save compiled data in a file.
%


function [data] = dyn_cell_packager(cellid, varargin)

p = inputParser();
addParameter(p, 'repack', false);
addParameter(p, 'datadir', []);
addParameter(p, 'dyn_path', []);
addParameter(p, 'do_save', true);
addParameter(p, 'file_prefix',  'cell_packager_data_');
addParameter(p, 'krn_width',   0.1); % kernel std (in secs)  
addParameter(p, 'krn_type', 'halfgauss'); % 'halfgauss' or 'fullgauss' for causal or acausal filter
addParameter(p, 'bin_size',  0.025); % (in secs)
parse(p, varargin{:});
par     = p.Results;

repack      = par.repack;
datadir     = par.datadir; 
do_save     = par.do_save;
file_prefix = par.file_prefix;
krn_width   = par.krn_width;
krn_type    = par.krn_type;
bin_size    = par.bin_size;

if isempty(par.datadir)
    if isempty(par.dyn_path)
        dp      = set_dyn_path;
    else
        dp      = par.dyn_path;
    end
    datadir = dp.celldat_dir;
end
    
filename = sprintf('%s%i.mat', file_prefix, cellid);
save_path = fullfile(datadir, filename);

if ~exist(datadir, 'dir'), mkdir(datadir); end;
if ~repack && exist(save_path, 'file')
	load(save_path);
	return;
end;

% get packaged data from cell. This includes event times and spike times
% without any further processing.
[array_data, vec_data, sessid, ratname] = package_dyn_phys(cellid,'save_to_file',0,'repack',repack);
[array_data] = compute_gen_state(array_data);
[array_data] = compute_state_switches(array_data);

data.cellid       = cellid;
data.ratname      = ratname;
data.sessid       = sessid;
% boolean that specifies in the recording was on right side of brain; if
% false, it was on left
data.brainsideright = bdata('select recorded_on_right from cells where cellid="{S}"', cellid);
[data.ntot_trials, data.sessiondate] = ...
	bdata('select n_done_trials, sessiondate from sessions where sessid="{S}"', sessid);
data.sessiondate  = data.sessiondate{1};
data.ngood_trials = numel(vec_data.good);
% smoothing and binning info
data.krn_width  = krn_width;
data.bin_size   = bin_size;
data.krn_type   = krn_type;

eibid = bdata('select eibid from cells where cellid="{S}"',cellid);

% assigns a region if it is from the following "regions" list.
% NOTE: For this code to work, the string listed in the EIBs table must
% contain a case-insensitive version of one of the items from the following
% list and it must not contain more than item from that list. 
regions = {'ppc','fof','sc','striatum'};
for i=1:numel(regions)
    region_eibs = bdata('select eibid from ratinfo.eibs where region like "{S}"',['%' regions{i} '%']);
    if(ismember(eibid, region_eibs))
        data.region = regions{i};
    end
end


% Calculate individual trial rate functions and get spike counts for different alignments

% Aligments are specified in the function called on the next line: align_LUT
    [align_strs, align_args] = dyn_align_LUT(2);
% loop through each alignment and get smoothed rate functions and spike counts
for i=1:numel(align_strs)
    [data.frate{i}, data.frate_t{i}] = make_rate_functions(cellid, 'array_data', array_data, 'vec_data', vec_data,...
                                                            'krn_width',krn_width,'bin_size',bin_size,align_args{i}{:},...
                                                            'krn_type', krn_type);
    [data.sp_counts{i}, data.sp_count_T{i}] = calc_sp_counts(cellid, 'array_data', array_data, 'vec_data', vec_data, align_args{i}{:});
end
data.align_strs = align_strs;
data.align_args = align_args;

% calculate normalization value
% it is mean firing rate across all trials at time = 0 aligned to stim_start
stimstart_index = strmatch('stimstart',align_strs,'exact');
[~,norm_time_bin] = min(abs(data.frate_t{stimstart_index}));
data.norm_mean = nanmean(data.frate{stimstart_index}(:,norm_time_bin));
data.norm_std = nanstd(data.frate{stimstart_index}(:,norm_time_bin));  

% % calculate normalized firing rates for each alignment
% % NOTE: not adding this because it doesn't seem the space-speed tradeoff makes it worth it
% for i=1:numel(align_strs)
%     data.fratenorm{i} = data.frate{i} ./ data.norm_mean;
%     data.fratez{i} = (data.frate{i} - data.norm_mean) ./ data.norm_std;
% end


ntrials = data.ngood_trials;
cpoke_start = vec_data.cpoke_start;

% data.trials is a sub-structure that keeps trial-by-trial data. It has
% fields that each have length equal to n-trials.
data.trials = struct('gamma', zeros(ntrials,1),...
    'rate_1', zeros(ntrials,1), 'rate_2', zeros(ntrials,1),...
    'lpulses', {cell(ntrials,1)}, ...
	'rpulses',     {cell(ntrials,1)},  ...
	'bup_diff', zeros(ntrials,1),  'bup_diff_rate',  zeros(ntrials,1), ...
	'spike_times', {cell(ntrials,1)},  'correct_dir',  zeros(ntrials,1), ...
	'rat_dir',     zeros(ntrials,1),   'hit',          zeros(ntrials,1), ...
	'stim_end',   zeros(ntrials,1), ...
    'cpoke_end',   zeros(ntrials,1),   'stim_start',   zeros(ntrials,1), ...
	'rt',          zeros(ntrials,1),   'T',            zeros(ntrials,1), ...
    'cpoke_out',   zeros(ntrials,1),   'hazard',       zeros(ntrials,1),...
    'genEndState', zeros(ntrials,1), ...
    'genSwitchTimes',{cell(ntrials,1)},  ...
    'switch_to_0', {cell(ntrials,1)}, 'switch_to_1', {cell(ntrials,1)});



for i=1:ntrials,
	data.trials.gamma(i)       = vec_data.gamma(i);
	data.trials.lpulses{i}     = array_data(i).left_bups(:)'  - cpoke_start(i);
	data.trials.rpulses{i}     = array_data(i).right_bups(:)' - cpoke_start(i);
    data.trials.bup_diff(i)    = vec_data.bup_diff(i);
	data.trials.spike_times{i} = array_data(i).spikes - cpoke_start(i);
	data.trials.correct_dir(i) = sign(data.trials.bup_diff(i));
	data.trials.rat_dir(i)     = 2*vec_data.pokedR(i) - 1; % +1 for R, -1 for L
	if data.trials.correct_dir(i) ~= 0,
		data.trials.hit(i)     = double(data.trials.correct_dir(i) == data.trials.rat_dir(i));
	else
		data.trials.hit(i)     = 0.5;   % assign 0.5 if there is no correct answer
	end;
	data.trials.cpoke_end(i)   = vec_data.cpoke_end(i)  - cpoke_start(i);
	data.trials.cpoke_out(i)   = vec_data.cpoke_out(i)  - cpoke_start(i);
	data.trials.stim_end(i)    = data.trials.cpoke_out(i);
    data.trials.stim_start(i)  = vec_data.stim_start(i) - cpoke_start(i);
	data.trials.rt(i)          = vec_data.spoke_in(i) - vec_data.cpoke_end(i);	
	data.trials.T(i)           = data.trials.cpoke_end(i) - data.trials.stim_start(i);
    data.trials.bup_diff_rate(i)  = data.trials.bup_diff(i) / data.trials.T(i);
    data.trials.rate_1(i) = array_data(i).rate_1;
    data.trials.rate_2(i) = array_data(i).rate_2;
    data.trials.hazard(i) =  array_data(i).hazard;
    data.trials.genEndState(i) = array_data(i).genEndState;
    data.trials.genSwitchTimes{i} = array_data(i).genSwitchTimes;
    data.trials.switch_to_0{i} = array_data(i).switch_to_0;
    data.trials.switch_to_1{i} = array_data(i).switch_to_1;

end;

% find preferred side and p-value of effect for each alignment
% it is choice side with higher mean "spike count" firing rate 
hitbool = data.trials.hit == 1;
for i=1:numel(align_strs)
    spct_rates = data.sp_counts{i} ./ data.sp_count_T{i};   % number of spikes divided by time
    right_mn = nanmean(spct_rates(vec_data.pokedR));        % mean of this for right choices
    left_mn = nanmean(spct_rates(~vec_data.pokedR));        % mean of this for left choices
    
    hit_mn = nanmean(spct_rates(hitbool));
    err_mn = nanmean(spct_rates(~hitbool));
    
    if right_mn > left_mn
        data.prefsideright{i} = 1;
    else
        data.prefsideright{i} = 0;
    end
    if hit_mn > err_mn
        data.prefoutcomehit{i} = 1;
    else
        data.prefoutcomehit{i} = 0;
    end
    % stats
    [~, data.prefp{i}] = ttest2(spct_rates(vec_data.pokedR), spct_rates(~vec_data.pokedR));
    [data.prefroc{i}, data.prefp_roc{i}] = bootroc(spct_rates(vec_data.pokedR), spct_rates(~vec_data.pokedR));
    
    [~, data.hitprefp{i}] = ttest2(spct_rates(hitbool), spct_rates(~hitbool));
    [data.hitprefroc{i}, data.hitprefp_roc{i}] = bootroc(spct_rates(hitbool), spct_rates(~hitbool));
end


% Calcualte time-dependent p-values
% This is done for pairs of alignments and condition sorting as specified
% in the function stats_LUT. 
[stats_strs, stats_aligns, stats_conds] = dyn_stats_LUT;
for i=1:numel(stats_strs)
    stat_ind = stats_aligns(i);   
    conds = eval(stats_conds{i});
    data.p_t{i} = data.frate_t{stat_ind};
    data.p{i} = calc_pvalue(data.frate{stat_ind}, conds, 'stats', 'ttest');
%     data.p{i} = calc_pvalue(data.frate{stat_ind}, conds, 'stats', 'roc');
    data.min_p{i} = min(data.p{i});
end
data.stats_strs = stats_strs;
data.stats_aligns = stats_aligns;
data.stats_conds = stats_conds;

% compute "lapse rate": performance on longest (top 25%) and easiest (highest gamma) trials
% Note that this is called "lapse rate", but it is really an estimate of
% the inverse of lapse rate using the performance on the super easy
% trials.
ag = unique(abs(data.trials.gamma));
sorted_T = sort(data.trials.T); 
topfourth = sorted_T(ceil(numel(sorted_T)*3/4));    % used to select longest trials
mytrials = find(abs(data.trials.gamma)==max(ag) & data.trials.T >= topfourth);
data.lapse_rate       = mean(data.trials.hit(mytrials));
data.lapse_rate_n     = numel(mytrials);
data.lapse_rate_t     = topfourth;  
data.lapse_rate_gamma = max(ag);    % gamma used to calculate "lapse rate"
data.min_gamma = min(ag);       % the smallest gamma included in data

% saves packaged data as mat file so that it doesn't need re-packaging
if do_save
    save(save_path, 'data');
end;


