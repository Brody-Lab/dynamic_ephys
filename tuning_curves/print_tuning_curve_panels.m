clear
close
% This script produces the figures from figure 3 or 5
% To produce panels from figure 3, set use_switches = 0
% To produce panels from figure 5, set use_switches = 1
use_switches = 0;
% which_switch determines whether to use 'model' switches or 'generative'
% switches. This is only relevant if use_switches = 1
which_switch = 'model';

if use_switches
    switch_params   = struct('t_buffers', [.2 .2]);
else
    which_switch    = '';
    switch_params   = [];
end

% set all the parameters that govern the tuning curves
lw          = 2;
fht         = 2.5;
fw          = 1.4*fht;
ppos        = [0 10 fw fht ];
fig_type    = '-dsvg';
do_print    = 1;
krn_width   = 0.1;
force_frdv  = 1;            % keep as one if rebinning
force_bin   = 0;
force_dv    = 1;
norm_type   = 'none';
dvlims      = [];
dp          = set_dyn_path;
direction   = 'backward';
krn_type    = 'halfgauss';
end_mask_s  = 0.0;
alignment   = 'stimstart';
dv_bins     = linspace(-6.5,6.5,9);
aticks      = [-7.5 -2.5 2.5 7.5];
lag         = 0.1;
max_dur     = 1;
switch_t0s  = [-.55:.025:.55];

switch_str      = '';
min_switch_t    = 0;
max_switch_t    = 3;
zscore_frates   = 0;
demean_frates   = 0;

do_normalize_tuning     = 1;

use_fake                    = 0;
do_shuffle_trials           = 0;
do_shuffle_fixing_choice    = 0;
do_shuffle_switches         = 0;
model_ratnoise              = [];
model                       = [];


ops.mask_after_stim_firing    = 1;
ops.mask_stim_firing          = 0;
ops.average_a_vals            = 1;
ops.average_fr                = 1;
ops.min_fr_mod                = 1;
ops.fit_time_sigmoid          = 0;
ops.use_nresidual             = 1;
ops.plot_res                  = 3;


%% analyze example cell 18181
excellid            = 18181;
data                = dyn_cell_packager(excellid);

switch alignment
    case 'stimend'
        max_dur     = 2;
        t0s         = -max_dur:0.025:-0.1-lag;  %
    case 'stimstart'
        t0s         = 0.1:0.025:max_dur-lag;  %
end

if ~isempty(which_switch)
    t0s_        = switch_t0s;
    switch_str  = ['_' which_switch(1:5) 'switch'];
    badtind = t0s_ > -.175 & t0s_ < .175;
    xlab = ['Time from ' which_switch ' switch (s)'];
else
    t0s_    = t0s;
    badtind = t0s_ > Inf;
    switch alignment
        case 'stimstart'
            xlab    = 'Time from stim on (s)';
        case 'stimend'
            xlab    = 'Time from stim off (s)';
    end
end

if ~use_fake
    res_fn = @(do_shuffle_switches, model, t0s_) ...
        dyn_fr_dv_map(excellid, 't0s', t0s_, ...
        'data',data, 'model', model, 'lag', lag, ...
        'alignment', alignment,...
        'krn_width', krn_width, 'krn_type', krn_type,...
        'norm_type', norm_type,'shuffle_trials',0,...
        'n_dv_bins',dv_bins,'end_mask_s',end_mask_s,...
        'demean_frates',demean_frates,'zscore_frates',zscore_frates,...
        'which_switch',which_switch, ...
        'switch_params',switch_params,...
        'shuffle_trials', do_shuffle_trials, ...
        'shuffle_switches', do_shuffle_switches,...
        'min_switch_t',min_switch_t, 'max_switch_t',max_switch_t,...
        'shuffle_trials_fixing_choice',do_shuffle_fixing_choice);
    print_fn = @(fh, str, fig_type) print(fh,fullfile(dp.fig_dir,[str switch_str]),...
        fig_type,'-painters');
else
    fname = fullfile(dp.model_dir,'simulate_model_ratnoise');
    load(fname, 'asample_ratnoise','model_ratnoise','t');
    fr_gen = -sign(asample_ratnoise);
    res_fn = @(do_shuffle_switches) dyn_fr_dv_map(18181, 't0s', t0s_, ...
        'frates', fr_gen, 'ft', t,...
        'model', model_ratnoise, 'data', data,...
        'lag', lag, 'alignment', alignment,...
        'krn_width', krn_width, 'krn_type', krn_type,...
        'norm_type', norm_type,'shuffle_trials',0,...
        'n_dv_bins',dv_bins,'end_mask_s',0,...
        'which_switch',which_switch,...
        'switch_params',switch_params,...
        'demean_frates',demean_frates,'zscore_frates',zscore_frates,...
        'do_svd',1, 'shuffle_trials', do_shuffle_trials, ...
        'shuffle_switches', do_shuffle_switches,...
        'min_switch_t',min_switch_t, 'max_switch_t',max_switch_t,...
        'shuffle_trials_fixing_choice',do_shuffle_fixing_choice);
    print_fn = @(fh, str, fig_type) print(fh,fullfile(dp.fig_dir,...
        ['fake_' str switch_str]),fig_type,'-painters');
end

res = res_fn(do_shuffle_switches, model, t0s_);
if ~do_print
    print_fn = @(fh, str, fig_type) fh;
end

fprintf('\nsvd rank 1 explains %.1f %% of the variance\n',100*res.rank_var(1))
%% Main tuning curve figure
if do_normalize_tuning
    ra = res.rank1_ra_n;
    ra = ra - min(ra);
    ra = ra ./ max(ra);
    res.rank1_ra_nn = ra;
end

if ~isempty(which_switch)
    title_str = sprintf('cell %i $E[r|a,t-t_c]$', res.cellid);
else
    title_str = sprintf('cell %i $E[r|a,t]$', res.cellid);
end

ra_field    = 'rank1_ra_nn';
mt_field    = 'rank1_mt_n';
ylab        = '\Delta FR';
mt          = [res.(mt_field)]';
ra          = [res.(ra_field)]';

% example cell E[r(a,t)]
fh = figure(1); clf
set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]));
ax = axes;
[~, cb] = fgta_line_plot(res, 'goodtind',~badtind, ...
    'linewidth', lw, 'ax', ax, ...
    'plot_field','fr_given_ta','dvlims',dvlims);
title(title_str,'fontweight','normal', 'interpreter', 'latex');
xlabel(xlab);
if demean_frates
    ylabel('\Delta FR (spikes/s)');
else
    ylabel('FR (spikes/s)');
end
xlim(res.t0s([1 end]));
axpos = get(ax,'position');
title(cb,sprintf('a'));
d = 1/length(res.dv_axis);
cb.Ticks = [0 1] + [1 -1]*(d*1.5);
cb.Position = cb.Position + [0.02 .15 -.02 -.4];
cb.TickLabels = round(10*res.dv_axis([2 end-1]))/10;

set(ax,'position',axpos);
print_fn(fh,'example_fgta',fig_type);
%% example residual
fh = figure(2); clf
ppos = [0 6 fw fht ];
set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]));
ax_r = axes;
fgta_line_plot(res,'goodtind',~badtind,...
    'ax',ax_r,'plot_field','fgta_resid','dvlims',dvlims);
xlabel(xlab)
ylabel({'\Delta FR (spikes/s)'})
if ~isempty(which_switch)
    title_str = sprintf('cell %i $E[\\Delta r|a,t-t_c]$', res.cellid);
    ylim(ax_r, [-1 1]*10)
else
    title_str = sprintf('cell %i $E[\\Delta r|a,t]$', res.cellid);
    ylim(ax_r, [-1 1]*7)
end
title(title_str,'fontweight','normal','interpreter','latex')

xlim(res.t0s([1 end]))
hold on
plot([0 0],ylim,'-k','linewidth',.1,'color',[1 1 1].*.5)
axpos = get(gca,'position');
print_fn(fh,'example_fgta_resid',fig_type);

%% print average tuning curve
fh = figure(4); clf
ppos = [0 2 fw fht ];
set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],...
    'papersize',ppos([3 4]));
axt = axes;
set(axt,'fontsize',12);
fgta_plot_tuning(res,'plot_field','fga_resid_tmn', ...
    'errorbar_field','fga_resid_std','ax',axt,'dvlims',dvlims,...
    'linewidth', 1);
pbaspect(axt, [1 1 1]);

xlabel('accumulated value (a)');
ylabel({'\Delta FR (spikes/s)'});
xlim(res.dv_axis([1 end])+[-1 1]);
set(axt,'XTick',aticks);
hold on
ylim([-5 5]+[-.5 .5]);
plot([0 0],ylim,'-k');
tfsz = ax_r.Title.FontSize;
title('$E[\Delta r|a]$','fontsize',tfsz,'fontweight','normal','interpreter','latex')
axt.Position([2 4]) = axpos([2 4]);
print_fn(fh,'example_fga_tmn',fig_type);

%% plot rank 1 approximation
ppos    = [5 6 fw fht ];
res     = compute_rank1_fgta_approx(res);
fh      = figure(3); clf
ax      = axes;

set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]));
fgta_line_plot(res,'goodtind',~badtind,'ax',ax,'plot_field','map_hat','dvlims',dvlims);
title('rank 1 approximation');
xlabel(xlab);
ylabel({'\Delta FR (spikes/s)'});
axis(ax, 'tight');
ylim(get(ax_r,'ylim'));
print_fn(fh,'example_rank1_map',fig_type);

res.rank1_mt_n2 = res.rank1_mt_n / 2;
res.rank1_ra_n2 = res.rank1_ra_n * 2;
ra_field    = 'rank1_ra_n';
mt_field    = 'rank1_mt_n';
ra_ylab     = {'normalized FR'};
mt_ylab     = 'FR modulation (spikes/s)';
%% Print rank 1 tuning curve
fh      = figure(5); clf
ppos    = [5 2 fw fht ];
ax      = axes;
set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]));

fgta_plot_tuning(res,'plot_field',ra_field, ...
    'errorbar_field','','ax',ax,'dvlims',dvlims, 'linewidth', 1);
pbaspect(ax, [1 1 1]);
title('rank 1 $\hat{f}(a)$','interpreter','latex');
xlabel('accumulated value (a)');
ylabel(ra_ylab);
box off
xlim(res.dv_axis([1 end])+[-1 1])
ylim([-.55 .55])
hold on
plot([0 0],ylim,'-k')
set(ax,'XTick',aticks,'TickDir', 'out')
print_fn(fh,'example_rank1_tuning',fig_type);
%% Print rank 1 modulation
ppos    = [10 2 fw fht ];
fh      = figure(6); clf
ax      = axes;
x       = res.t0s;
set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]))

if ~isempty(badtind)
    xx          = res.t0s;
    xx(badtind) = nan;
    yy          = res.(mt_field);
    yy(badtind) = nan;
    bad0        = find(badtind,1)-1;
    badn        = find(badtind,1,'last')+1;
    plot(xx, yy, 'k', 'linewidth', 1.5)
    hold on
    plot(xx([bad0 badn]), yy([bad0 badn]), ':k', 'linewidth', 1)
    plot([0 0],ylim,':k')
else
    plot(res.t0s, res.(mt_field), 'k', 'linewidth', 1)
    hold on
end
if isempty(which_switch)
    plot([.5 .5],ylim,':k')
end
xlim(x([1 end]))
if ~isempty(which_switch)
    title_str = 'rank 1 $\hat{m}(t-t_c)$';
else
    title_str = 'rank 1 $\hat{m}(t)$';
end
title(title_str,'interpreter','latex');
xlabel(xlab);
ylabel(mt_ylab);
box off;
cb = colorbar;
drawnow
pos = get(ax,'position');
delete(cb);
set(ax,'position',pos, 'TickDir', 'out');
print_fn(fh,'example_rank1_mod',fig_type);

%% population level tuning curves
% load up population level results
refit   = 0;
frbins  = [];
demean_frates = 0;
zscore_frates = 1;
switch which_switch
    case ''
        fun = @(pop_cellids,pop_sessids) pop_dyn_fr_dv_map(pop_cellids, pop_sessids,...
            't0s', t0s, 'lag', lag, ...
            'alignment', alignment,...
            'krn_width', krn_width, 'krn_type', krn_type,...
            'norm_type', norm_type,'shuffle_trials',0,...
            'demean_frates',demean_frates,'zscore_frates',zscore_frates,...
            'n_dv_bins',dv_bins,'end_mask_s',end_mask_s);
        
        fig_prefix = '';
    case 'model'
        fun = @(pop_cellids,pop_sessids) pop_dyn_fr_dv_map(pop_cellids, pop_sessids,...
            't0s', switch_t0s, 'lag', lag, ...
            'alignment', alignment,...
            'krn_width', krn_width, 'krn_type', krn_type,...
            'norm_type', norm_type,'shuffle_trials',0,...
            'n_dv_bins',dv_bins,'end_mask_s',end_mask_s,...
            'demean_frates',demean_frates,'zscore_frates',zscore_frates,...
            'which_switch','model',...
            'switch_params',switch_params);
        fig_prefix = 'modelswitch_';
    case 'generative'
        fun = @(pop_cellids,pop_sessids) pop_dyn_fr_dv_map(pop_cellids, pop_sessids,...
            't0s', switch_t0s, 'lag', lag, ...
            'alignment', alignment,...
            'krn_width', krn_width, 'krn_type', krn_type,...
            'norm_type', norm_type,'shuffle_trials',0,...
            'n_dv_bins',dv_bins,'end_mask_s',end_mask_s,...
            'demean_frates',demean_frates,'zscore_frates',zscore_frates,...
            'which_switch','generative',...
            'switch_params',switch_params);
        fig_prefix = 'genswitch_';
end

sp = get_switch_params(which_switch,switch_params);
if zscore_frates
    fig_prefix = [fig_prefix 'zscore_'];
end
if isempty(which_switch)
    fn = fullfile(dp.model_fits_dir, ['pop_' fig_prefix...
        'tuning_res.mat']);
else
    stadir = get_sta_dirname(sp);
    fn = fullfile(stadir, ['pop_' fig_prefix...
        'tuning_res.mat']);
end
if refit | ~exist(fn,'file')
    % get the list of good cells
    cout_auc_file   = fullfile(dp.ephys_summary_dir,'cout_auc.mat');
    v               = load(cout_auc_file)
    pop_cellids     = v.cellids(v.good_cells);
    pop_cellids     = pop_cellids(1:end);
    pop_sessids     = v.sessids;
    pop_res         = fun(pop_cellids,pop_sessids);
    save(fn, 'pop_res','-v7.3')
else
    load(fn,'pop_res')
end

%% flip cell maps by preference
do_normalize_tuning = 0; % don't need to normalize this time because we already zscored
ra_field    = 'rank1_ra_n';
mt_field	= 'rank1_mt_n';
mt          = [pop_res.(mt_field)]';
ra          = [pop_res.(ra_field)]';

if do_normalize_tuning
    ra      = ra - min(ra,[],2);
    ra      = ra ./ max(ra,[],2);
    ylab    = {'fractional change' 'in firing rate r'}
else
    ylab    = '\Delta FR';
end

dv_axis     = pop_res(1).dv_axis;
flip_cell   = mean(ra(:,dv_axis>0),2) < mean(ra(:,dv_axis<0),2);

pop_fgta    = cat(3,pop_res.fr_given_ta);
pop_fgta_r  = cat(3,pop_res.fgta_resid);

ra(flip_cell,:) = fliplr(ra(flip_cell,:));

pop_fgta_r(:,:,flip_cell)   = fliplr(pop_fgta_r(:,:,flip_cell));
pop_fgta(:,:,flip_cell)     = fliplr(pop_fgta(:,:,flip_cell));

mean_pop_fgta_r = mean(pop_fgta_r,3);
mean_pop_fgta   = mean(pop_fgta,3);

tempres         = pop_res(1);
tempres.fgta_resid  = mean_pop_fgta_r;
tempres.fr_given_ta = mean_pop_fgta;

tempres = compute_rank1_fgta_approx(tempres,'which_map','fr_given_ta');
tempres.dv_axis = dv_axis;

%% plot variance explained histogram
all_var = horzcat(pop_res.rank_var);
all_var = all_var(1,:);
fh = figure(3); clf
set(fh,'position',[5 5 fht fht])
histogram(all_var)
fprintf(['\nfor single cells, svd rank 1 explains an average of '...
    '%.1f %% of the variance SD %.1f%% \n'],...
    100*mean(all_var),100*std(all_var))

fprintf('population average svd rank 1 explains %.1f %% of the variance\n',...
    100*tempres.rank_var(1))

%% Plot population average map
pref_clr = [.8 .25 .8];
npref_clr = [.8 .65 .25];
fh = figure(21); clf

set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]))
ax = axes;
fgta_line_plot(tempres,'goodtind',~badtind,'ax',ax,...
    'plot_field','fr_given_ta','dvlims',dvlims,...
    'up_clr', pref_clr, 'down_clr', npref_clr);
ylabel(ax,'\Delta FR (z-score)')
xlabel(ax,xlab)
title('population average','interpreter','latex')
if zscore_frates
    ylim([-1 1].*.13)
    set(ax, 'ytick', [-.1 0 .1])
end
print_fn(fh,'pop_fgta',fig_type)
%% population average residual map
fh      = figure(20); clf
ppos    = [1 1 fw fht ];
set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]))
ax      = axes;

[ax, cb] = fgta_line_plot(tempres,'goodtind',~badtind,'ax',ax,...
    'plot_field','fgta_resid','dvlims',dvlims,...
    'up_clr', [.8 .25 .8], 'down_clr', [.8 .65 .25]);
ylabel(ax,'\Delta FR (z-score)')
title('average of selective cells','interpreter','latex')
xlabel(xlab)
if zscore_frates
    if isempty(which_switch)
        ylim([-1 1].*.125)
    else
        ylim([-1 1].*.15)
    end
    set(ax,'ytick',[-.1 0 .1]);
else
    ylim([-1 1].*max(abs(ylim)))
end
title(cb,{'a for ' 'pref. side'})
d = 1/length(tempres.dv_axis)
drawnow
axpos = get(ax,'position')
cb.Ticks = [0 1] + [1 -1]*(d*1.5);
cb.Position = cb.Position + [0.02 .15 -.02 -.4];
cb.TickLabels = round(10*tempres.dv_axis([2 end-1]))/10;

set(ax,'position',get(ax_r,'position'))
print_fn(fh,'pop_fgta_resid',fig_type)

%% POP TUNING 
this_tmn = tempres.(ra_field);

if do_normalize_tuning
    this_tmn    = tempres.(ra_field);
    this_tmn_   = this_tmn - min(this_tmn)
    tempres.(ra_field) = this_tmn_;
end
mt_ylab = 'FR modulation (z-score)';

fh = figure(7); clf
ppos = [5 2 fw fht ];
set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],...
    'papersize',ppos([3 4]))
ax = axes;
fgta_plot_tuning(tempres,'linecolor','k',...
    'plot_field',ra_field, ...
    'errorbar_field','','ax',ax,'dvlims',dvlims,'linewidth',1,...
    'up_clr',pref_clr, 'dwn_clr', npref_clr);
hold(ax,'on')
pbaspect(ax, [1 1 1])
set(ax,'XTick',aticks)
title('rank 1 $\hat{f}(a)$','interpreter','latex')
xlabel('accumulated value (a)')
ylabel(ra_ylab)
if zscore_frates
    ylim([-.55 .55])
end
%
plot([0 0],ylim,'-','color',[1 1 1]*.8)
box off
ax.TickDir = 'out';
drawnow
xlim(res.dv_axis([1 end])+[-1 1])
xx = tempres.dv_axis;
this_b = fit_four_param_psycho(xx(:),this_tmn(:));
slope_pop_mean = prod(this_b([2 3]))/4;
print_fn(fh,'pop_ra',fig_type)
%% POP MODULATION
fh = figure(8); clf
ppos = [5 5 fw fht ];
set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],...
    'papersize',ppos([3 4]));
ax = axes;
hold on
if ~isempty(badtind)
    xx = tempres.t0s;
    xx(badtind) = nan;
    yy = tempres.(mt_field);
    yy(badtind) = nan;
    plot(xx, yy, 'k', 'linewidth', 1.5)
    bad0 = find(badtind,1)-1;
    badn = find(badtind,1,'last')+1;
    hold on
    plot(xx([bad0 badn]), yy([bad0 badn]), ':k', 'linewidth', 1)
    ylim([0 max(ylim)])
    plot([0 0],ylim,'-', 'color', [1 1 1].*.8)
else
    plot(tempres.t0s, tempres.(mt_field), 'k', 'linewidth', 1)
    hold on
    plot([.5 .5],ylim, '-k', 'linewidth', 1)
end
title('rank 1 $\hat{m}(t)$','interpreter','latex')
xlabel(xlab)
box off
ylim([0 max(ylim)])
xlim(x([1 end]))
ylabel(mt_ylab)
pos = get(ax_r,'position');
set(ax,'position',pos, 'TickDir', 'out');

print_fn(fh,'pop_frm',fig_type)





%% plot all selective tuning curves and modulations - not used
%close all

fignamefun  = @(x) fullfile(dp.fig_dir,[fig_prefix x]);
print_fn    = @(fh,x,fig_type) print(fh,fignamefun(x),fig_type,'-painters');

fh = figure(1); clf
set(fh,'position',[5 5 fw fht])
set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]))
x = pop_res.t0s;
plot(x([1 end]),[0 0],'k')
hold on
plot(x,mt,'color',[1 1 1].*.7)
plot(x,mean(mt),'k','linewidth',2)
box off

xlabel(xlab)
xlim(x([1 end]))
cb = colorbar
if zscore_frates
    ylabel({'modulation' '(z-scores)'})
end

drawnow
pos = get(gca,'position')

delete(cb)
set(gca,'position',pos)
print_fn(fh,['selective_' mt_field],fig_type)
 

fh = figure(11); clf
set(fh,'position',[5 5 fht fht])
xx = pop_res.dv_axis;
plot(xx,ra,'color',[.5 .5 .5])
hold on
plot(xx,mean(ra),'k','linewidth',2)

box off
ylabel(ylab)
xlabel('accumulated value (a)')
xlim(tempres.dv_axis([1 end])+[-1 1])
if ~do_normalize_tuning & strcmp(ra_field, 'rank1_ra_n')
    ylim([-.61 .61])
elseif do_normalize_tuning
    ylim([-.1 1.1])
end
%ylim([-5.5 5.5])
set(gca, 'TickDir', 'out')
print_fn(fh,['selective_' ra_field],fig_type)
%% plot modulation as heat map - not used
fh = figure(11); clf
ppos = [5 5 fw fht ];
set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]))
ax = axes
[~, sortind] = sort(mean(mt,2),'descend');
imagesc(mt(sortind,:),'x',pop_res(1).t0s)
caxis([0 .3])
colormap(bone)
cb = colorbar('northoutside');
title(cb,'modulation (z-scores)')
xlabel('Time from model switch')
ylabel('cell #')
print_fn(fh,'pop_frm_heat',fig_type)

%% plot tuning curves as heat map - not used
mt_n = [pop_res.rank1_mt_n]';
ra_n = [pop_res.rank1_ra_n]';

ra_n(flip_cell,:) = fliplr(ra_n(flip_cell,:));
mean_mt = mean(mt_n,2);
mean_ra_n = mean_mt .* ra_n;
[~, sortorder] = sort(mean_mt,'descend');
figure(111); clf
imagesc(mean_ra_n(sortorder,:),'x',dv_axis)
colormap(colormapRedBlue)
caxis([-2 2]*.1)
colorbar('location','southoutside')
%% make a histogram of slopes - not used
%pop_fga_tmn = vertcat(pop_res.fga_tmn);
fitfun = @(beta, clickdiff) fourParamPsychometric(beta, clickdiff);
p0 = [.1 .8 .01 0]; % [left lapse (right lapse - left lapse) slope controller bias]
opts.MaxIter = 100000;
opts.TolFun = 1e-15;
opts.TolX = 1e-15;
xx = pop_res(1).dv_axis;
slope_field = 'rank1_ra_n';
tmns = horzcat(pop_res.(slope_field))';
figure(1); clf
for ii = 1:length(pop_res)
    %%
    this_tmn = pop_res(ii).(slope_field) ;
    if mean(this_tmn(xx<0)) > mean(this_tmn(xx>0))
        this_tmn = -this_tmn;
    end
    this_tmn = this_tmn - min(this_tmn);
    this_tmn = this_tmn ./ max(this_tmn);
    tmns(ii,:) = this_tmn';
    
    b(ii,:) = fit_four_param_psycho(xx(:),this_tmn(:))
    
    subplot(121); cla
    plot(xx,fourParamPsychometric(b(ii,:), xx(:)),'r');
    hold on
    plot(xx,this_tmn,'k.-');
    ylim([0 1])
    
    subplot(122);
    
    hold on
    plot(xx,this_tmn,'.-','color',[.5 .5 .5]);
    ylim([0 1])
    drawnow
    
    %%
    %pause()
end
plot(xx,mean(tmns),'k','linewidth',2)

slopes = b(:,3).*b(:,2)/4;
fh = figure(10); clf
ppos = [1 6 fw fht ]
set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]))
histogram(slopes,[0:.025:.325],'facecolor',colors('purple'))
ylabel('number of neurons')
xlabel('slope')
box off
hold on
plot([1 1].*slope_pop_mean,ylim,'k')
%plot([1 1].*mean(slopes),ylim,'color',colors('purple'))
plot(percentile(slopes,[5 95]),[1 1]*max(ylim),'-','color',colors('purple'))
plot(mean(slopes),[1 1]*max(ylim),'.','color',colors('purple'))
print_fn(fh,['pop_slopes' slope_field],fig_type)
%%

