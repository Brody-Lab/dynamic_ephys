function [real_diff, diff_p, mt_p] = post_switch_gain(varargin)
p = inputParser;
addParameter(p, 'which_switch', 'model');
addParameter(p, 'alpha', .05);
addParameter(p, 'n_shuffles', 250);
addParameter(p, 't_buffers', [0 0]);
addParameter(p, 'model_smooth_wdw', 100);
parse(p,varargin{:});
p = p.Results;

switch_params   = struct('which_switch',p.which_switch,...
        't_buffers',t.buffers,...
        'model_smooth_wdw',p.model_smooth_wdw,...
        'clear_bad_strengths', 1, 'bad_strength', 0, 'fit_line', 1,...
        'min_pre_dur', 0, 'min_post_dur', 0,...
        'min_switch_t',0,'max_switch_t',Inf,...
        'exclude_final',0,'final_only',0);
    
dp              = set_dyn_path;
cout_auc_file   = fullfile(dp.ephys_summary_dir,'cout_auc.mat');
v               = load(cout_auc_file,'cellids','good_cells','sessids');
pop_cellids     = v.cellids(v.good_cells);
pop_sessids     = v.sessids(v.good_cells);
unique_sessid   = unique(pop_sessids);

this_stadir     = get_sta_dirname(switch_params);
if ~exist(this_stadir,'dir')
    mkdir(this_stadir);
end
fn = fullfile(this_stadir,['gain_change.mat']);
%%
alpha       = p.alpha;
nc          = length(pop_cellids);
recompute   = 0;
fprintf('working through %i sessions', length(unique_sessid));

%% RERUN ALL CELLS
n_shuffles  = p.n_shuffles;
real_diff   = nan(nc,1);
diff_p      = nan(nc,1);
counter     = 0;
for si = 1:length(unique_sessid)
    fprintf('session %i...',si)
    this_sessid = unique_sessid(si);
    ind         = find(pop_sessids == unique_sessid(si));
    for ii = 1:length(ind)
        cc          = ind(ii);
        counter     = counter + 1;
        excellid    = pop_cellids(cc);
        fprintf('\ncell %i...',counter)
        data        = dyn_cell_packager(excellid);
        if ii == 1
            model   = 	get_data_model_p(data, data.trials.trialnums);
        end
        
        [real_diff(cc), diff_p(cc), mt_p_] = switch_gain_shuffle(cellid, n_shuffles,...
            'data',data,'model',model,'switch_params',switch_params);
        if counter == 1
            mt_p = nan(nc,size(mt_p_,2));
        end
        mt_p(cc,:) = mt_p_;
    end
end

sig_hi  = diff_p > (1-alpha/2);
sig_lo  = diff_p < alpha/2;

save(fn,'real_diff','diff_p','sig_hi','sig_lo')

