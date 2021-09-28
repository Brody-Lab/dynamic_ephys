cellid = 18181;
which_switch = 'model';

res =  compute_switch_triggered_average(cellid, 'force', 1,...
    'which_switch', which_switch, 'n_shuffles', 0, ...
    'save_file', 0,'min_pre_dur', .4, 'min_post_dur', .4,...
    'include_str', 'true(size(data.trials.hit == 1))',...
    't_buffers',[.2 .2]);
STR_left    = res.STR_left_real;
STR_right   = res.STR_right_real;
lags        = res.lags;

%plotSTA(res,'lag',.1)
%%
STR_left    = res.STR_left_real;
STR_right   = res.STR_right_real;
max_t = .55;
min_t = -.55;
cx = [-1 1].*1*max([STR_left(:); STR_right(:)]);
cx = [-1 1].* 50;
%cx = [20 80]

dp = set_dyn_path;

fht = 2.5;
fw = 4.375; %1.75*2.5;
fht = 2.5;
fw  = 3.75;
fh = figure(1); clf
set(fh,'position',[2 2 fw fht],'papersize', [fw fht])
axh = axes;

n = 15;
off = 0;%-200;
%indA = randi(size(STR_left,1),1,n);% + [0:-1:-n];
%indB = randi(size(STR_right,1),1,n);% + [0:-1:-n];
indA = -off + size(STR_left,1) + [0:-1:-n+1];
indB = -off + size(STR_right,1) + [0:-1:-n+1];

A = STR_left(indA,:);
B = STR_right(indB,:);

NA = ~isnan(A);
NB = ~isnan(B);

ax1 = subplot(211);
ax2 = subplot(212);

bg_color = [1 1 1] .* .85;

imagesc((A),'x',lags,'Alphadata',NA,'parent',ax1)
title(ax1,'1 $\rightarrow$ 2','interpreter','latex')
title(ax1,'go right $\rightarrow$ go left','interpreter','latex')

imagesc((B),'x',lags,'Alphadata',NB,'parent',ax2)
title(ax2,'2 $\rightarrow$ 1','interpreter','latex')
title(ax2,'go left $\rightarrow$ go right','interpreter','latex')

set(ax1,'color',bg_color)
set(ax2,'color',bg_color)

set([ax1 ax2], 'xlim',[min_t max_t])
colormap(colormapRedBlue)
cb = colorbar;
drawnow 
cb.Position = cb.Position + [.125 .325 -.02 -.125];
%ax2pos = get(ax2,'position');
%set(ax2,'position',ax2pos)

title(cb,'\Delta fr')
ax1.TickDir = 'out';
box(ax1,'off')
set(ax1,'ytick',[], 'xtick', [])
ax2.TickDir = 'out';
box(ax2,'off')
set(ax2,'ytick',[])
xlabel(ax2,['Time from ' which_switch ' state change (s)']) 
ylabel(ax2,'state switch #') 
ylabel(ax1,'state switch #') 
hold(ax1,'on')
hold(ax2,'on')
plot(ax1,[0 0], ylim, 'k--')
plot(ax2,[0 0], ylim, 'k--')
caxis(ax1,cx)
caxis(ax2,cx)

print(fh, fullfile(dp.fig_dir, 'example_sta_trials'), '-dsvg','-painters')
