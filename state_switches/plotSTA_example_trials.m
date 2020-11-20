cellid = 18181;
which_switch = 'model';
norm_type = 'z';
res = plotSTA(cellid,which_switch,1);
%%
res =  compute_switch_triggered_average(cellid,...
    'post', 2, 'force', 1, 'final_only',0,...
   'which_switch',which_switch, 'n_shuffles', 100,...
   'save_file',1,'mask_other_switch',0, 'compute_residual',1);
STR_left    = res.STR_left_real;
STR_right   = res.STR_right_real;
lags        = res.lags;
%%
% % 
% % ind = 1:30;
dp = set_dyn_path;
  max_t = .75;
    min_t = -.45;
fh = figure(1); clf
set(fh,'position',[2 2 6 3],'papersize', [6 3])
axh = axes;
% max_y = max([STR_left(:); STR_right(:)]);
% for ii = 1:3
%     subplot(121); hold on
%     plot(lags,ii+STR_left(ii,:)./max_y,'k');
%     subplot(122); hold on
%     plot(lags,ii+STR_right(ii,:)./max_y,'k');
% end

n = 15;
indA = size(STR_left,1) + [0:-1:-n];
indB = size(STR_right,1) + [0:-1:-n];

% indA = 1:100;
% indB = 1:100;

A = (STR_left(indA,:));
B = (STR_right(indB,:));

NA = ~isnan(A);
NB = ~isnan(B);

ax1 = subplot(211)
ax2 = subplot(212)

bg_color = [1 1 1] .* .75;

n = 60
indA = randi(size(STR_left,1),n,1);
indB = randi(size(STR_left,1),n,1);

imagesc((A),'x',lags,'Alphadata',NA,'parent',ax1)
title(ax1,'2 $\rightarrow$ 1','interpreter','latex')
imagesc((B),'x',lags,'Alphadata',NB,'parent',ax2)
title(ax2,'1 $\rightarrow$ 2','interpreter','latex')


set(ax1,'color',bg_color)
set(ax2,'color',bg_color)

cx = [-1 1].*50
caxis(ax1,cx)
caxis(ax2,cx)

set([ax1 ax2], 'xlim',[min_t max_t])
colormap(colormapRedBlue)
drawnow 
ax2pos = get(ax2,'position');
cb = colorbar;
cb.Position = cb.Position + [.11 0 -.01 -.15];
set(ax2,'position',ax2pos)
title(cb,'\Delta fr')
ax1.TickDir = 'out';
box(ax1,'off')
set(ax1,'ytick',[], 'xtick', [])
ax2.TickDir = 'out';
box(ax2,'off')
set(ax2,'ytick',[])
xlabel(ax2,'time from state switch (s)') 
ylabel(ax2,'# state switch') 
hold(ax1,'on')
hold(ax2,'on')
plot(ax1,[0 0], ylim, 'k--')
plot(ax2,[0 0], ylim, 'k--')

% colormap(dyn_cmap(1000))
% cm = colormapLinear(dp.left_color);
% cm = [flipud(cm); colormapLinear(dp.right_color)];
% colormap(cm)
print(fh, fullfile(dp.fig_dir, 'example_sta_trials'), '-dsvg','-painters')
%%
ax3 = subplot(322)
ax4 = subplot(324)
ax5 = subplot(325)

plot(ax5,lags,res.dprime_shuff,'color',[1 1 1].*.9)
hold(ax5,'on')
plot(ax5,lags,res.dprime_real,'r')

hold(ax3,'on')
hold(ax4,'on')
plot(ax3, lags, (STR_left(indA,:)),'color',bg_color,'linewidth',.5)
plot(ax3, lags, nanmean(A),'color',dp.left_color,'linewidth',.5)

plot(ax4, lags, (STR_right(indB,:)),'color',bg_color,'linewidth',.5)
plot(ax4, lags, nanmean(B),'color',dp.right_color,'linewidth',.5)
%plot(ax3, lags, (B),'color',dp.right_color,'linewidth',.5)
%%
plotSTA(18181,'model',0)
