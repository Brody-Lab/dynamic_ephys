cellid = 18181;
data = dyn_cell_packager(cellid);
model = get_data_model_p(data);
%%

% having a hard time understanding why using only hits or only errors but
% not both would remove the dip in the middle

lag = 0.1
dt = .025
hits_only = 0;
errs_only = 1;
which_trials = [];
pre_dur = .1;
post_dur = .1;
switch 1
    case 0
        t0s = .1:dt:1.25;
        which_switch = '';
    case 1
        t0s = -.5:dt:.5;
        which_switch = 'model';
end
dv_bins = [-6:.5:6];

res = dyn_fr_dv_map(cellid, 'lag', lag, ...
    't0s',t0s,'n_dv_bins',dv_bins,'which_switch',which_switch,...
    'model',model,'data',data,...
    'hits_only',hits_only,'errs_only',errs_only,...
    'which_trials',which_trials,...
    'min_pre_dur', pre_dur, 'min_post_dur', post_dur);
%%
dmn_res = dyn_fr_dv_map(cellid, 'lag', lag, ...
    'n_dv_bins',dv_bins,'demean_frates',1,'do_svd',0,...
    't0s',t0s ,'which_switch',which_switch,...
    'model',model,'data',data,...
    'hits_only',hits_only,'errs_only',errs_only,...
    'whicH_trials',which_trials,...
    'min_pre_dur', pre_dur, 'min_post_dur', post_dur);

z_res = dyn_fr_dv_map(cellid, 'lag', lag, ...
    'n_dv_bins',dv_bins,'norm_type','div','do_svd',0,...
    't0s',t0s ,'which_switch',which_switch,...
    'model',model,'data',data,...
    'hits_only',hits_only,'errs_only',errs_only,...
    'whicH_trials',which_trials,...
    'min_pre_dur', pre_dur, 'min_post_dur', post_dur);
%%
mvmnsz = 1;

fh = figure(1); clf
set(fh,'position',[1 7 4 8])
subplot(211)
fgta_line_plot(res,'plot_field','fr_given_ta','mvmnsz',mvmnsz)
subplot(212)
fgta_line_plot(res,'mvmnsz',mvmnsz)

fh = figure(2); clf
set(fh,'position',[6 7 4 8])
subplot(211)
fgta_line_plot(dmn_res,'plot_field','fr_given_ta','mvmnsz',mvmnsz)
subplot(212)
fgta_line_plot(dmn_res,'mvmnsz',mvmnsz)

fh = figure(3); clf
set(fh,'position',[11 7 4 8])
subplot(211)
fgta_line_plot(z_res,'plot_field','fr_given_ta','mvmnsz',mvmnsz)
subplot(212)
fgta_line_plot(z_res,'mvmnsz',mvmnsz)
%%
figure(4); clf
res.mass_ta = squeeze(sum(res.mass_tfa,2))./res.mass_t/10;

use_val = res.mass_ta > .01;

subplot(211)
imagesc(res.mass_ta','x',res.t0s,'y',res.dv_axis)
colormap(gca,flipud(bone))
colorbar
subplot(212)
fgta_line_plot(res,'plot_field','mass_ta')
hold on
%plot(res.t0s,res.mass_t,'k')


%%
res.fr_given_ta(~use_val) = nan;
res.fgta_resid(~use_val) = nan;
mvmnsz = 1;
ptile = 0;
dvlims = [-5 5]/2
fh = figure(1); clf
set(fh,'position',[1 7 4 8])
subplot(211)
fgta_line_plot(res,'plot_field','fr_given_ta','mvmnsz',mvmnsz,...
    'plot_mass_ptile',ptile,'dvlims',dvlims)
subplot(212)
fgta_line_plot(res,'mvmnsz',mvmnsz,...
    'plot_mass_ptile',ptile,'dvlims',dvlims)


fh = figure(2); clf
set(fh,'position',[6 7 4 8])
subplot(211)
fgta_line_plot(dmn_res,'plot_field','fr_given_ta','mvmnsz',mvmnsz,...
    'plot_mass_ptile',ptile,'dvlims',dvlims)
subplot(212)
fgta_line_plot(dmn_res,'mvmnsz',mvmnsz,...
    'plot_mass_ptile',ptile,'dvlims',dvlims)

fh = figure(3); clf
set(fh,'position',[11 7 4 8])
subplot(211)
fgta_line_plot(z_res,'plot_field','fr_given_ta','mvmnsz',mvmnsz,...
    'plot_mass_ptile',ptile,'dvlims',dvlims)
subplot(212)
fgta_line_plot(z_res,'mvmnsz',mvmnsz,...
    'plot_mass_ptile',ptile,'dvlims',dvlims)