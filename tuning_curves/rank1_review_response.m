dp = set_dyn_path(1);
fht = 2.5;
fw = 1.75*fht;
ppos = [28 6 fw fht ];
lw = 2;
dv_bins = linspace(-6.5,6.5,9);
dvlims  = dv_bins([1 end]);
%%
% load tuning curves for included population
%fig_prefix = 'modelswitch_';
which_switch    = '';

fig_prefix = '';
fig_prefix = [fig_prefix 'zscore_'];
ra_ylab = '\Delta FR (z-score)';
fn = fullfile(dp.model_fits_dir, ['pop_' fig_prefix...
    'tuning_res.mat']);

load(fn,'pop_res')
%% plot VE for ranks 1-5

fh = figure(1); clf
rank_var = [pop_res(:).rank_var]';
figure(1); clf 
subplot(121)
histogram(rank_var(:,1),20,'facecolor',[1 1 1].*.5)
hold on
plot([1 1].*mean(rank_var(:,1)),ylim,'k')

ylabel('# cells')
xlabel('variance explained by rank 1')
box off
subplot(122)
plot(1:size(rank_var,2), rank_var, '-','color',[1 1 1].*.8,'linewidth',1)
hold on
errorbar(1:size(rank_var,2), mean(rank_var), 1.96*sem(rank_var),...
    'ko-','capsize',0,'linewidth',2,'markerfacecolor','w')
box off
ylabel('variance explained')
xlabel('rank')
%% plot worst cell
pop_cellids = [pop_res(:).cellid];
[low_ve, low_ve_ind] = min(rank_var(:,1));
[sorted_ve, ve_sort_ind] = sort(rank_var(:,1));
nc = length(sorted_ve);
ii = 3;
low_ve_ind = ve_sort_ind(ii);

res     = pop_res(low_ve_ind);
t0s_    = res.t0s;
alignment   = 'stimstart';

if ~isempty(which_switch)
    badtind = t0s_ > -.175 & t0s_ < .175;
    xlab = ['time from ' which_switch ' switch (s)'];
    title_str = sprintf('cell %i $E[r|a,t-t_c]$', res.cellid);
else
    badtind = t0s_ > Inf;
    switch alignment
        case 'stimstart'
            xlab  = 'time from stim on (s)';
        case 'stimend'
            xlab = 'time from stim off (s)';
    end
    title_str = sprintf('cell %i $E[r|a,t]$', res.cellid);
end


fh = figure(2); clf
set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]));

ax = axes;
[~, cb] = fgta_line_plot(res,'goodtind',~badtind, ...
    'linewidth', lw, 'ax', ax, ...
    'plot_field','fr_given_ta','dvlims',dvlims);

fh = figure(3); clf
set(fh,'position',ppos+[5 0 0 0],'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]));

ax_r = axes;
fgta_line_plot(res,'goodtind',~badtind,...
    'ax',ax_r,'plot_field','fgta_resid','dvlims',dvlims);

fh = figure(4); clf
set(fh,'position',ppos+[0 -4 0 0],'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]));

ax = axes;
fgta_line_plot(res,'goodtind',~badtind,'ax',ax,'plot_field','map_hat','dvlims',dvlims)

fh = figure(5); clf
set(fh,'position',ppos+[5 -4 0 0],'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]));
ax = axes;
res = compute_rank1_fgta_approx(res);
rn = 5;
res.map_hat_5 = (res.u(:,1:rn)'.*res.s(1:rn))'*res.v(:,1:rn)';
fgta_line_plot(res,'goodtind',~badtind,'ax',ax,...
    'plot_field','map_hat_5','dvlims',dvlims)
title(sprintf('rank %i approximation',rn))

ra_field = 'rank1_ra_n';
mt_field = 'rank1_mt_n';

fh = figure(6); clf

set(fh,'position',ppos + [0 4 0 0],'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]))
ax = axes;
fgta_plot_tuning(res,'plot_field',ra_field, ...
    'errorbar_field','','ax',ax,'dvlims',dvlims, 'linewidth', 1)
pbaspect(ax, [1 1 1])
title('rank 1 $\hat{f}(a)$','interpreter','latex')
xlabel('accumulated value (a)')

ylabel(ra_ylab)
box off
ax.TickDir = 'out';
xlim(res.dv_axis([1 end])+[-1 1])
ylim([-.55 .55])
hold on
plot([0 0],ylim,'-k')
set(ax,'XTick',aticks)

mt_ylab = 'FR modulation (z-score)';
fh = figure(7); clf
set(fh,'position',ppos + [5 4 0 0],'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]))
ax = axes;
x = res.t0s;
if ~isempty(badtind)
    xx = res.t0s;
    xx(badtind) = nan;
    yy = res.(mt_field);
    yy(badtind) = nan;
    plot(xx, yy, 'k', 'linewidth', 1.5)
    bad0 = find(badtind,1)-1;
    badn = find(badtind,1,'last')+1;
    hold on
    plot(xx([bad0 badn]), yy([bad0 badn]), ':k', 'linewidth', 1)
    plot([0 0],ylim,':k')
else
    plot(res.t0s, res.(mt_field), 'k', 'linewidth', 1)
    hold on
end
if ~isempty(which_switch)
    plot([.5 .5],ylim,':k')
end
xlim(x([1 end]))

if ~isempty(which_switch)
title_str = 'rank 1 $\hat{m}(t-t_c)$';
else
    title_str = 'rank 1 $\hat{m}(t)$';
end
box off
plot(xlim,[0 0],'k:','linewidth',1)
    title(title_str,'interpreter','latex')
xlabel(xlab)
ylabel(mt_ylab)
box off
ax.TickDir = 'out';
cb = colorbar
drawnow 
pos = get(ax,'position')
delete(cb)
set(ax,'position',pos)
%print_fn(fh,'example_rank1_mod',fig_type);

aticks      = [-7.5 -2.5 2.5 7.5];
fh = figure(8); clf
set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],...
    'papersize',ppos([3 4]));
axt = axes;
set(axt,'fontsize',12)
fgta_plot_tuning(res,'plot_field','fga_resid_tmn', ...
    'errorbar_field','fga_resid_std','ax',axt,'dvlims',dvlims,...
    'linewidth', 1)
pbaspect(axt, [1 1 1])

xlabel('accumulated value (a)')
ylabel({'\Delta FR (spikes/s)'})
xlim(res.dv_axis([1 end])+[-1 1]);
set(axt,'XTick',aticks);
hold on
%ylim([-5 5]+[-.5 .5]);
plot([0 0],ylim,'-k');
%tfsz = ax_r.Title.FontSize;
title('$E[\Delta r|a]$','fontsize',14,'fontweight','normal','interpreter','latex')
%ylim(ax, [-1 1]*7)
%axt.Position([2 4]) = axpos([2 4]);
% print_fn(fh,'example_fga_tmn',fig_type);
% plot a psth
example_cell_psth('cells',res.cellid, 'fpos', ppos + [-4.5 0 0 0])
%%
