do_print        = 1;
dv_bins         = linspace(-6.5,6.5,9);

aticks      = [-7.5 -2.5 2.5 7.5];
lag         = 0.1;            %0.2;
max_dur     = 1;
switch_t0s  = [-.55:.025:.55];
n_dv_bins   = 100;          % triggers rebinning, but not recompiling 
krn_width   = 0.01;
force_frdv  = 1;            % keep as one if rebinning
force_bin   = 0;
force_dv    = 1;
norm_type   = 'none';
ops.mask_after_stim_firing    = 1;
ops.mask_stim_firing          = 0;
ops.average_a_vals            = 1;
ops.average_fr                = 1;
ops.min_fr_mod                = 1;
ops.fit_time_sigmoid          = 0;
ops.use_nresidual             = 1;
ops.plot_res                  = 3;
dvlims = [];
t0s         = 0.1:0.025:max_dur-lag;  %

ratname = 'H066';
fht = 2.5;
fw = 1.75*fht;
dp          = set_dyn_path;
direction   = 'backward';

krn_type    = 'halfgauss';
end_mask_s  = 0.0;
alignment   = 'stimstart'; %'cpokeout'

excellid = 18362;
excellid = 18699;
data = dyn_cell_packager(excellid);
[~, vec_data, ~,~] = get_behavior_data(dp.spikes_dir, ...
    data.cellid, data.sessid,ops);
[model, constant_x, allsame] = get_data_model_p(data, vec_data);
[align_strs, align_args] = dyn_align_LUT;
align_ind = strmatch(alignment,align_strs,'exact');
nanfrates   = isnan(data.frate{align_ind});
        
min_switch_t = 0;
max_switch_t = 3;
fig_type = '-dsvg';
norm_type = '';
use_fake = 0;
zscore_frates = 0;
demean_frates = 0;
switch_str = '';
do_shuffle_trials = 0;
do_shuffle_fixing_choice = 0;
do_shuffle_switches = 0;
%%
nt = length(data.trials.trialnums);
xv1 = false(nt,1);
ind = randperm(nt);
xv1(ind(1:ceil(nt/2))) = true;

xv2 = ~xv1;


which_switch = '';
res_xv1 = dyn_fr_dv_map(excellid, 'which_trials', xv1,...
    't0s', t0s, ...
    'data',data, 'model', model, 'lag', lag, ...
    'alignment', alignment,...
    'krn_width', krn_width, 'krn_type', krn_type,...
    'norm_type', norm_type,'shuffle_trials',0,...
    'n_dv_bins',dv_bins,'end_mask_s',end_mask_s,...
    'demean_frates',demean_frates,'zscore_frates',zscore_frates,...
    'which_switch',which_switch, 'shuffle_trials', do_shuffle_trials, ...
    'shuffle_switches', do_shuffle_switches,...
    'min_switch_t',min_switch_t, 'max_switch_t',max_switch_t,...
    'shuffle_trials_fixing_choice',do_shuffle_fixing_choice);

res_xv2 = dyn_fr_dv_map(excellid, 'which_trials', xv2,...
    't0s', t0s, ...
    'data',data, 'model', model, 'lag', lag, ...
    'alignment', alignment,...
    'krn_width', krn_width, 'krn_type', krn_type,...
    'norm_type', norm_type,'shuffle_trials',0,...
    'n_dv_bins',dv_bins,'end_mask_s',end_mask_s,...
    'demean_frates',demean_frates,'zscore_frates',zscore_frates,...
    'which_switch',which_switch, 'shuffle_trials', do_shuffle_trials, ...
    'shuffle_switches', do_shuffle_switches,...
    'min_switch_t',min_switch_t, 'max_switch_t',max_switch_t,...
    'shuffle_trials_fixing_choice',do_shuffle_fixing_choice);

figure(1); clf
ax1 = subplot(211)
fgta_line_plot(res_xv1, 'ax', ax1)

ax2 = subplot(212)
fgta_line_plot(res_xv2, 'ax', ax2)

linkaxes([ax1 ax2],'y')
