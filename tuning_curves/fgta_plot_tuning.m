function [ax, plotres] = fgta_plot_tuning(res,varargin)
p = inputParser;
addParameter(p,'linewidth',1)
addParameter(p,'clrs',[])
addParameter(p,'dwn_clr',[])
addParameter(p,'up_clr',[])
addParameter(p,'ax',[])
addParameter(p,'mvmnsz',1)
addParameter(p,'plot_field','fga_tmn')
addParameter(p,'errorbar_field','fga_sem')
addParameter(p,'plot_colorbar',1)
addParameter(p,'dvlims',[])
addParameter(p,'linecolor','k')
addParameter(p,'linestyle','-')
addParameter(p,'markersize',6)
parse(p, varargin{:})
p = p.Results;


dp = set_dyn_path;

if isempty(p.ax)
    ax = axes;
else
    ax = p.ax;
end

if isempty(p.dvlims)
    p.dvlims = [res.dv_axis(1) res.dv_axis(end)];
end

lw          = p.linewidth;
mvmnsz      = p.mvmnsz;
plot_mean   = res.(p.plot_field);
dv_keep     = res.dv_axis >= p.dvlims(1) & res.dv_axis <= p.dvlims(end);
dv_axis     = res.dv_axis(dv_keep);
if all(dv_keep([1 end]) == 1) & length(unique(diff(dv_axis)))==2
    special_endpoints = 1;
else
    special_endpoints = 0;
end
n_dv_bins   = length(dv_axis);
cellid      = res.cellid;
t0s         = res.t0s;

if ~isempty(p.dvlims)
    plot_mean = plot_mean(dv_keep);
    if ~isempty(p.errorbar_field)
        plot_errbar  = res.(p.errorbar_field);
        plot_errbar  = plot_errbar(dv_keep);
    end
end

if isempty(p.clrs)
    mid_color = [1 1 1].*.875;
    nc = floor((n_dv_bins)/2);
    if isempty(p.up_clr)
        p.up_clr = dp.right_color;
    end
    if isempty(p.dwn_clr)
        p.dwn_clr = dp.left_color;
    end
    cm = flipud(colormapLinear(p.dwn_clr,nc,mid_color));
    cm2 = colormapLinear(p.up_clr,nc,mid_color);
    if mod(n_dv_bins,2) == 0
    clrs = [cm(2:end,:); cm2(1:end-1,:)];
    else 
    clrs = [cm(1:end-1,:); cm2(1:end,:)];
    end
    assert(size(clrs,1)==n_dv_bins);
else
    clrs = p.clrs;
end


if mvmnsz == 1
    smooth_fn = @(x) x;
    smooth_clrs = clrs;
else
    smooth_fn   = @(x) movmean(movmean(x,mvmnsz),mvmnsz,2);
    smooth_clrs = clrs;
    % smooth_fn   = @(x) movmean(x(:,1:mvmnsz:end),mvmnsz);
    % smooth_clrs = clrs(:,1:mvmnsz:end);
end
hold(ax, 'on')
if special_endpoints
    plot(dv_axis([2:end-1]),plot_mean([2:end-1]),p.linestyle,'color',p.linecolor,'linewidth',lw,...
    'markersize',p.markersize)

plot(dv_axis([1 2]),plot_mean([1 2]),':','color',p.linecolor,'linewidth',lw,...
    'markersize',p.markersize)
plot(dv_axis([end-1 end]),plot_mean([end-1 end]),':','color',p.linecolor,'linewidth',lw,...
    'markersize',p.markersize)
else
plot(dv_axis,plot_mean,p.linestyle,'color',p.linecolor,'linewidth',lw,...
    'markersize',p.markersize)
end

plotres.dv_axis = dv_axis;
plotres.plot_mean = plot_mean;

for i=1:mvmnsz:numel(dv_axis)
    hold on;
    if ~isempty(p.errorbar_field)
        eh = errorbar(ax,dv_axis(i), plot_mean(i), ...
            plot_errbar(i), '.', 'markersize',10,'Color', clrs(i,:),'linewidth',2,...
            'capsize',0);
        plotres.plot_errbar(i) = plot_errbar(i);
    else
        eh = plot(ax,dv_axis(i), plot_mean(i), ...
        '.', 'markersize',10,'Color', clrs(i,:),'linewidth',2);
    end;
end
%ylim(max(abs(plot_mean))*[-1 1])
%fr_mod = 0.75;
%plot([results.dv_axis(1) results.dv_axis(end)], [fr_mod/2 fr_mod/2], 'r--');
%plot([results.dv_axis(1) results.dv_axis(end)], [-fr_mod/2 -fr_mod/2], 'r--');
title('raw r(a)')
xlabel('accumulated value')
ylabel('normalized FR')
ax.TickDir = 'out';
box(ax,'off')
% plot([dv_axis(1) dv_axis(end)], [0 0], 'k--');
% 
% ylims = ylim;
% ylim([-max(abs(ylims)), +max(abs(ylims))])
% plot([0 0], [-max(abs(ylims)), +max(abs(ylims))], 'k--');