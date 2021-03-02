function ax = fgta_heat_plot(res,varargin)
p = inputParser;
addParameter(p,'linewidth',2)
addParameter(p,'clrs',[])
addParameter(p,'ax',[])
addParameter(p,'mvmnsz',1)
addParameter(p,'plot_field','fgta_resid')
addParameter(p,'plot_colorbar',1)
addParameter(p,'dvlims',[])
addParameter(p,'plot_mass_ptile',[])
parse(p, varargin{:})
p = p.Results;

dp = set_dyn_path;

if isempty(p.ax)
    ax = gca;
else
    ax = p.ax;
end

if isempty(p.dvlims)
    p.dvlims = [res.dv_axis(1) res.dv_axis(end)];
end

lw          = p.linewidth;
mvmnsz      = p.mvmnsz;
to_plot     = res.(p.plot_field);
dv_keep     = res.dv_axis >= p.dvlims(1) & res.dv_axis <= p.dvlims(end)
dv_axis     = res.dv_axis(dv_keep);
n_dv_bins   = length(dv_axis);
n_fr_bins   = length(res.frbins);
cellid      = res.cellid;
t0s         = res.t0s;

mass_ta     = squeeze(sum(res.mass_tfa,2))./res.mass_t/10;
if ~isempty(p.plot_mass_ptile) & p.plot_mass_ptile > 0
    doplot      = mass_ta > percentile(mass_ta,p.plot_mass_ptile);
else
    doplot = true(size(mass_ta));
end



if ~isempty(p.dvlims)
    to_plot = to_plot(:,dv_keep);
    doplot  = doplot(:,dv_keep);
end

to_plot(~doplot) = nan;

% if isempty(p.clrs)
%     mid_color = [1 1 1].*.925;
%     nc = floor((n_fr_bins)/2);
%     cm = flipud(colormapLinear(dp.left_color,nc,mid_color));
%     cm2 = colormapLinear(dp.right_color,nc,mid_color);
%     if mod(n_fr_bins,2) == 0
%     clrs = [cm(2:end,:); cm2(1:end-1,:)];
%     else 
%     clrs = [cm(1:end-1,:); cm2(1:end,:)];
%     end
%     assert(size(clrs,1)==n_fr_bins);
% else
%     clrs = p.clrs;
% end

clrs = colormapRedBlue;

if mvmnsz == 1
    smooth_fn = @(x) x;
    smooth_clrs = clrs;
else
    smooth_fn   = @(x) movmean(movmean(x,mvmnsz),mvmnsz,2);
    smooth_clrs = clrs;
    % smooth_fn   = @(x) movmean(x(:,1:mvmnsz:end),mvmnsz);
    % smooth_clrs = clrs(:,1:mvmnsz:end);
end

% set(ax,'ColorOrder',smooth_clrs,'NextPlot','ReplaceChildren',...
%     'fontsize',12)
% 
% set(ax,'ColorOrder',smooth_clrs,'NextPlot','ReplaceChildren',...
%     'fontsize',12)
%plot(ax,t0s,smooth_fn(to_plot),'-','linewidth',lw)
imagesc(ax,smooth_fn(to_plot)','x',t0s,'y',dv_axis)
ylabel('accumulated value (a)')
xlabel('time (s)')
title(['Cell ' num2str(cellid) ' r(a,t)'])

colormap(ax,clrs)
if p.plot_colorbar
    cb = colorbar;
end

set(ax,'TickDir','out')
axis xy
xlim(ax, res.t0s([1 end])+[-.1 .1])
