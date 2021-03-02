function axh = plot_cell_map(res, varargin)
p = inputParser;
addParameter(p, 'mvmnsz', 1)
addParameter(p, 'fh', [])
addParameter(p, 'linecolor', 'k')
parse(p,varargin{:})
p = p.Results;

dp = set_dyn_path;

mvmnsz  = p.mvmnsz;
if isempty(p.fh)
    fh = figure;
elseif isnumeric(p.fh)
    fh = figure(p.fh);
elseif ishandle(p.fh)
    fh = p.fh;
end



savefig = 0;
plot_sigmoid = 1;
plot_nice = 1;

markersize = 2;

doflip = 0;
lw = 2;

fga_tmn_n   = res.fga_tmn_n;
fgta        = res.fr_given_ta;
fgta_resid  = res.fgta_resid;
fgta_residn = res.fgta_resid_n;
fga_rnta    = res.fga_rn_tmn;
fga_std_nta = res.fga_rn_std;
if isfield(res,'map_hat')
rank1_map   = res.map_hat;
rank1_fa    = res.rank1_ra_n;
rank1_mt    = res.rank1_mt_n;
end
t0s         = res.t0s;
cellid      = res.cellid;
if isfield(res, 'svd_betas')
    svd_betas   = res.svd_betas;
else
    svd_betas = nan;
end
dv_x        = res.dv_axis;
fga_std_n   = res.fga_std_n;
fga_tmn     = res.fga_tmn;
fga_std     = res.fga_std;
frm_time    = res.frm_time;
fga_aa_cell = nanmean(res.fr_given_ta,2);
fga_resid_tmn = res.fga_resid_tmn;
fga_resid_std = res.fga_resid_std;

nbetas      = nan;
betas       = nan;

% plot fga_ta_cell, with sigmoid for both cases
n_dv_bins   = numel(fga_tmn_n);
%clrs        = color_set(n_dv_bins);
mid_color = [1 1 1].*.925;
nc = floor((n_dv_bins)/2);
cm = flipud(colormapLinear(dp.left_color,nc,mid_color));
cm2 = colormapLinear(dp.right_color,nc,mid_color);
if mod(n_dv_bins,2) == 0
clrs = [cm(2:end,:); cm2(1:end-1,:)];
else 
clrs = [cm(1:end-1,:); cm2(1:end,:)];
end
assert(size(clrs,1)==n_dv_bins);


if mvmnsz == 1
    smooth_fn = @(x) x;
    smooth_clrs = clrs;
else
    smooth_fn   = @(x) movmean(movmean(x,mvmnsz),mvmnsz,2);
    smooth_clrs = clrs;
    % smooth_fn   = @(x) movmean(x(:,1:mvmnsz:end),mvmnsz);
    % smooth_clrs = clrs(:,1:mvmnsz:end);
end


%%%% Plot r(a,t)
ax = subplot(4,4,1); 
hold on; 
set(ax,'ColorOrder',smooth_clrs,'NextPlot','ReplaceChildren',...
    'fontsize',12)
plot(ax,t0s,smooth_fn(fgta),'-','linewidth',lw)
ylabel('Firing Rate (Hz)')
xlabel('Time (s)')
title(['Cell ' num2str(cellid) ' r(a,t)'])
%ylim([ inf])
axis('tight')

ax = subplot(4,4,5);
imagesc(fgta,'y',t0s,'x',dv_x)
colormap(ax,colormapRedBlue)
cb = colorbar('northoutside')
title(cb, 'Firing Rate (Hz)')
ylabel('Time (s)')
xlabel('Accumulated Value (a)')
axis('tight')

%%%% Plot Delta r(a,t)
ax = subplot(4,4,2);
set(ax,'ColorOrder',smooth_clrs,'NextPlot','ReplaceChildren',...
    'fontsize',12)
plot(ax,t0s,smooth_fn(fgta_resid),'-','linewidth',lw);
ylabel('\Delta FR (Hz)')
xlabel('Time (s)')
title('\Delta r(a,t)')
ylims = ylim;
ylim([-max(abs(ylims)), +max(abs(ylims))])
axis('tight')

ax = subplot(4,4,6);
imagesc(fgta_resid,'y',t0s,'x',dv_x)
colormap(ax,colormapRedBlue)
cb = colorbar('northoutside')
title(cb, '\Delta FR (Hz)')
ylabel('Time (s)')
xlabel('Accumulated Value (a)')
axis('tight')

%%%% Plot Rank 1 Delta r(a,t)
ax = subplot(4,4,3); cla(ax)
set(ax,'ColorOrder',smooth_clrs,'NextPlot','ReplaceChildren',...
    'fontsize',12)
plot(t0s, smooth_fn(rank1_map),'-','linewidth',lw);
ylabel('\Delta FR (Hz)')
xlabel('Time (s)')
title('\Delta r(a,t)-r1')
title('rank 1 \Delta r(a,t)')
axis('tight')

ax = subplot(4,4,7);
imagesc(rank1_map,'y',t0s,'x',dv_x)
colormap(ax,colormapRedBlue)
cb = colorbar('northoutside')
title(cb, 'rank 1 \Delta r(a,t)')
ylabel('Time (s)')
xlabel('Accumulated Value (a)')
axis('tight')

%%%% Plot Normalized Delta r(a,t)
ax = subplot(4,4,4);cla; hold on; set(gca, 'fontsize',12);
% if plot_nice
%     fgta_residn = movmean(fgta_residn,mvmnsz);
%     fgta_residn = fgta_residn - min(fgta_residn,[],2);
%     fgta_residn = fgta_residn./max(fgta_residn,[],2);
% end
set(ax,'ColorOrder',smooth_clrs,'NextPlot','ReplaceChildren',...
    'fontsize',12)
plot(ax,t0s,smooth_fn(fgta_residn), ...
    '-', 'MarkerSize', markersize,'linewidth',lw);
ylabel('Norm. FR')
xlabel('Time (s)')
title('\Delta r(a,t) norm.')
colormap(ax,clrs)
drawnow
pos = get(ax,'position')
cb = colorbar('eastoutside')
set(ax,'position',pos)
set(cb, 'ytick', [0 1], 'yticklabel', dv_x([1 end]))
title(cb,'a')

ax = subplot(4,4,8);
imagesc(fgta_residn,'y',t0s,'x',dv_x)
colormap(ax,colormapRedBlue)
cb = colorbar('northoutside')
title(cb, '\Delta r(a,t) norm.')
ylabel('Time (s)')
xlabel('Accumulated Value (a)')
axis('tight')

%%%% Plot tuning curve from avg method?
subplot(4,4,9); 
hold on; set(gca, 'fontsize',12);
plot(dv_x,fga_tmn,'k','linewidth',2,'color',p.linecolor)
h = gca;
for i=1:mvmnsz:numel(dv_x)
    hold on;
    eh = errorplot2(h, dv_x(i), fga_tmn(i), ...
        fga_std(i), 'Marker', '', 'Color', clrs(i,:));
end;
xfine = dv_x(1):0.1:dv_x(end);
% if ~isnan(betas)
%     ypred = dyn_sig(betas,xfine);
%     if doflip
%         ypred = flipdim(ypred,2);
%     end
%     if plot_sigmoid
%         plot(xfine, ypred, 'm','linewidth',2)
%     end
% end
title('Avg. r(a)')
xlabel('Accumulated Value')
ylabel('Norm. FR')
%ylim([-0.1 1.1])
plot([0 0], ylim, 'k--');
%plot([dv_x(1) dv_x(end)], [0 0], 'k--');
which_t = t0s > .25 & t0s < 1;
%plot(dv_x,mean(fgta(which_t,:)),'m')
xlim(dv_x([1 end]))

%%%% Plot un-normalized tuning curve
h = subplot(4,4,10);
hold on; set(gca, 'fontsize',12);
plot(dv_x,fga_resid_tmn,'k','linewidth',2,'color',p.linecolor)
for i=1:mvmnsz:numel(dv_x)
    hold on;
    eh = errorplot2(h,dv_x(i), fga_resid_tmn(i), ...
        fga_resid_std(i), 'Marker', '', 'Color', clrs(i,:));
end;
%fr_mod = 0.75;
%plot([results.dv_axis(1) results.dv_axis(end)], [fr_mod/2 fr_mod/2], 'r--');
%plot([results.dv_axis(1) results.dv_axis(end)], [-fr_mod/2 -fr_mod/2], 'r--');
plot([dv_x(1) dv_x(end)], [0 0], 'k--');
title('raw r(a)')
xlabel('Accumulated Value')
ylabel('\Delta FR (Hz)')
ylims = ylim;
ylim([-max(abs(ylims)), +max(abs(ylims))])
plot([0 0], [-max(abs(ylims)), +max(abs(ylims))], 'k--');
xlim(dv_x([1 end]))

%%%% Plot SVD tuning curve
subplot(4,4,11); 
hold on; set(gca, 'fontsize',12)
xfine = dv_x(1):0.1:dv_x(end);
if ~isnan(svd_betas)
    ypred = dyn_sig(svd_betas,xfine);
    if doflip
        ypred = flipdim(ypred,2);
    end
    if plot_sigmoid
        plot(xfine, ypred+min(rank1_fa), 'm','linewidth',2)
    end
end

plot(dv_x,rank1_fa,'k','linewidth',2,'color',p.linecolor)
title('SVD r(a)')
xlabel('Accumulation Value (a)')
ylabel('FR')
ylim([-0.75 .75])
set(gca, 'Ytick', [-.5 0 0.5], 'yticklabel',{'-.5' '0','0.5'})
plot([dv_x(1) dv_x(end)], [0 0], 'k--');
plot([0 0], ylim, 'k--');
xlim(dv_x([1 end]))

%%%% Plot normalized tuning curve
subplot(4,4,12);
hold on; set(gca, 'fontsize',12)
plot(dv_x,fga_rnta,'k','linewidth',2,'color',p.linecolor)
h = gca;
for i=1:mvmnsz:numel(dv_x)
    hold on;
    eh = errorplot2(h,dv_x(i), fga_rnta(i), fga_std_nta(i), 'Marker', '', 'Color', clrs(i,:));
end;
xfine = dv_x(1):0.1:dv_x(end);
if ~isnan(nbetas)
    ypred = dyn_sig(nbetas,xfine);
    if doflip
        ypred = flipdim(ypred,2);
    end
    if plot_sigmoid
        plot(xfine, ypred, 'm','linewidth',2)
    end
end
title('Norm. r(a)')
xlabel('Accumulation Value (a)')
ylabel('Norm. FR')
plot([dv_x(1) dv_x(end)], .5+[0 0], 'k--');
plot([0 0], [0 1], 'k--');
ylim([-0.1 1.1])
plot([0 0], [0 1], 'k--');
set(h, 'Ytick', [0 0.5 1], 'yticklabel',{'0','0.5', '1'})
xlim(dv_x([1 end]))

%%%% Plot residual fr
subplot(4,4,13); 
hold on;
plot(t0s, fga_aa_cell, 'k','linewidth',2,'color',p.linecolor)
plot(t0s, zeros(size(t0s)),'Color',[.75 .75 .75],'linewidth',1);
title('Avg. FR')
ylabel('Avg. FR.(Hz)')
xlabel('Time (s)');
set(gca, 'fontsize',12);

%%%% Plot firing rate modulation
subplot(4,4,14); hold on;
plot(t0s, frm_time, 'k','linewidth',2,'color',p.linecolor)
plot(t0s, zeros(size(t0s)),'Color',[.75 .75 .75],'linewidth',1);
title('FR. Mod')
ylabel('FR. Mod. (Hz)')
xlabel('Time (s)');
set(gca, 'fontsize',12);

%%%% Plot firing rate modulation
subplot(4,4,15); hold on;
plot(t0s, rank1_mt, 'k','linewidth',2,'color',p.linecolor)
plot(t0s, zeros(size(t0s)),'Color',[.75 .75 .75],'linewidth',1);
title('SVD FR. Mod')
ylabel('FR mod. (Hz)')
xlabel('Time (s)');
set(gca, 'fontsize',12);

subplot(4,4,16); hold on;
%plot(results.t0s, var_expl.*100, 'k','linewidth',2)
nrank = length(res.rank_var);
plot(1:nrank, 100*res.rank_var, 'k','linewidth',2,'color',p.linecolor)
%plot([1 nrank], [0 0],'Color',[.75 .75 .75],'linewidth',1);
ylabel('VE %')
xlabel('Rank');
set(gca, 'fontsize',12);
%axis('tight')
ylim(100*[.6 1])

if savefig
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 10 6];
    fig.PaperPositionMode = 'Manual';
    fig.PaperSize = [10 6];
    print(gcf,  ['../../figures/tuning/' num2str(cellid) ],'-dsvg')
end

axh = get(gcf,'children');




