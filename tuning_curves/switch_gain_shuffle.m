function [diff_p, mt_p, res] = switch_gain_shuffle(cellid, n_shuffles, varargin)
p = inputParser;
addParameter(p,'data',[])
addParameter(p,'model',[])
addParameter(p,'do_save',1)
addParameter(p,'recompute',0)
parse(p,varargin{:})
p = p.Results;

if ~isempty(p.data)
    data = p.data;
else
    data = dyn_cell_packager(cellid);
end
if ~isempty(p.model)
    model   = p.model;
else
    model   = get_data_model_p(data, data.trials.trialnums);
end

res = [];

ra_field    = 'rank1_ra_n';
mt_field    = 'rank1_mt_n';

alpha       = .05;
krn_width   = 0.1;
dv_bins     = linspace(-6.5,6.5,9);
switch_t0s  = [-.55:.025:.55]; %#ok<NBRAK>
badtind     = switch_t0s > -.2 & switch_t0s < .20;
bad0        = find(badtind,1)-1;
badn        = find(badtind,1,'last')+1;

lag             = .1;
which_switch    = 'model';
switch_params   = struct('which_switch',which_switch,...
            't_buffers',[.2 .2],...
            'clear_bad_strengths', 1, 'bad_strength', 0, 'fit_line', 1,...
            'min_pre_dur', 0, 'min_post_dur', 0,...
            'min_switch_t',0,'max_switch_t',Inf,...
            'exclude_final',0,'final_only',0,'model_smooth_wdw',100);

this_stadir     = get_sta_dirname(switch_params);
if ~exist(this_stadir,'dir')
    mkdir(this_stadir);
end
fn = fullfile(this_stadir,[cellid '.mat']);

if exist(fn,'file') & ~p.recompute
    load(fn)
    return
end

res_fn = @(cellid, data, model, do_shuffle_switches) ...
    dyn_fr_dv_map(cellid, 'data', data, 'model', model, ...
    'shuffle_switches', do_shuffle_switches,...
    't0s', switch_t0s,  'lag', lag, ...
    'n_dv_bins', dv_bins, ...
    'which_switch', which_switch, ...
    'switch_params',switch_params);
% compute for real data
res                 = res_fn(cellid, data, model, 0);
mt_real             = res.(mt_field);
mt_real(badtind)    = nan;
real_diff           = nanmean(mt_real(badn:end)) - ...
    nanmean(mt_real(1:bad0));
% compute shuffles
mt_shuffle  = nan(n_shuffles, numel(res.(mt_field)));
ra_shuffle  = nan(n_shuffles, numel(res.(ra_field)));

fprintf('\nbeginning shuffle...\n')  
parfor ss = 1:n_shuffles        % compute the unshuffled results
    if mod(ss,25) == 0
        fprintf('%i (of %i)...',ss,n_shuffles)
    end
    res_shuff = res_fn(cellid, data, model, 1);
    ra_shuffle(ss,:)    = res_shuff.(ra_field);
    temp                = res_shuff.(mt_field);
    temp(badtind)       = nan;
    mt_shuffle(ss,:)    = temp; 
end

shuff_diff  = nanmean(mt_shuffle(:,badn:end),2) - ...
    nanmean(mt_shuffle(:,1:bad0),2);
mt_p        = mean(mt_real' > mt_shuffle,1);
diff_p      = mean(real_diff > shuff_diff);
%%
if 1            
    save(fn,'mt_p','diff_p')
end
