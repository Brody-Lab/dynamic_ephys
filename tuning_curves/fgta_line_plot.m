function [ax cb] = fgta_line_plot(res,varargin)
p = inputParser;
addParameter(p,'linewidth',2)
addParameter(p,'clrs',[])
addParameter(p,'up_clr',[])
addParameter(p,'down_clr',[])
addParameter(p,'ax',[])
addParameter(p,'mvmnsz',1)
addParameter(p,'plot_field','fgta_resid')
addParameter(p,'colorbar_location','eastoutside')
addParameter(p,'dvlims',[])
addParameter(p,'goodtind',[])
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
dv_keep     = res.dv_axis >= p.dvlims(1) & res.dv_axis <= p.dvlims(end);
dv_axis     = res.dv_axis(dv_keep);
n_dv_bins   = length(dv_axis);
t0s         = res.t0s;
if length(res) > 1
    cellid = [];
    cellidstr = 'population';
else
    cellid      = res.cellid;
    cellidstr = ['cell ' num2str(cellid)];
end


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

if ~all(p.goodtind)
    badtind = ~p.goodtind;  
    t0s(badtind) = nan;
    to_plot(badtind,:) = nan;
    badx0 = res.t0s(find(badtind,1)) - diff(res.t0s([1 2]))/2;
    badxn =  res.t0s(find(badtind,1,'last')) + diff(res.t0s([1 2]))/2;
    
end




if isempty(p.clrs)
    mid_color = [1 1 1].*.875;
    nc = floor((n_dv_bins)/2);
    if ~isempty(p.up_clr)
        up_clr = p.up_clr;
    else
        up_clr = dp.right_color;
    end
    if ~isempty(p.down_clr)
        down_clr = p.down_clr;
    else
        down_clr = dp.left_color;
    end
    cm = flipud(colormapLinear(down_clr,nc,mid_color));
    cm2 = colormapLinear(up_clr,nc,mid_color);
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

set(ax,'ColorOrder',smooth_clrs,'NextPlot','ReplaceChildren',...
    'fontsize',12)

set(ax,'ColorOrder',smooth_clrs,'NextPlot','ReplaceChildren',...
    'fontsize',12)
plot(ax,t0s,smooth_fn(to_plot),'-','linewidth',lw)

ylabel('Firing Rate (Hz)')
xlabel('Time (s)')
title([cellidstr ' r(a,t)'])

if ~isempty(p.colorbar_location)
    
    colormap(smooth_clrs);
    drawnow
    pos = get(ax,'position');
    cb = colorbar(p.colorbar_location);
    set(cb, 'ytick', [0 1], 'yticklabel', round(dv_axis([1 end])*10)/10)
    title(cb,'a')
    %set(ax,'position',pos)
end



set(ax,'TickDir','out')

xlim(ax, res.t0s([1 end])+[-.1 .1])

axis(ax, 'tight')

if exist('badx0','var')
    draw_patch = 0;
    if draw_patch
        patch([badx0  badx0 badxn badxn],...
            .98*[min(ylim) max(ylim)  max(ylim) min(ylim)],...
            [1 1 1].*.95,'edgecolor','w')
    end
    this_ind = [find(badtind,1)-1 find(badtind,1,'last')+1];
    hold on
    pp = plot(ax,t0s(this_ind),smooth_fn(to_plot(this_ind,:)),...
        ':','linewidth',1);
    

end

