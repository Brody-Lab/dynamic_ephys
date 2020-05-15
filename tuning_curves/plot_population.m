function [] = plot_population(results,varargin)

savefig = 0;
if nargin >1
savefig = varargin{1};
disp('Saving figures')
end
res = 3;

%%%% plot first figure
colorfig = figure(1); clf;
h = gca;
set(h, 'fontsize',16)
hold on;

n_dv_bins   = numel(results.dv_axis);
n_cells     = size(results.fga_aa_cell,2);
clrs        = color_set(n_dv_bins);
Markersize  = 2;

for i=1:res:n_dv_bins
    plot(h,results.t0s,movmean(results.fga_population(:,i),3)-0.5,'-','Color',clrs(i,:),'Markerfacecolor',clrs(i,:),'MarkerSize',Markersize);
end
leg_loc                 = 'Northwest';
legend_str{1}           = 'Negative bound';
legend_str{n_dv_bins}   = 'Positive bound';
for i=2:n_dv_bins-1
    legend_str{i} = ['DV = ' num2str(results.dv_axis(i))];
end
xlabel(h,'Time (s)');
ylabel(h,'Mean normalized FR');
title(h,'Population Average');
%ylim([0 1])
ylim([-.5 .5])
cpos = get(colorfig, 'Position');
pbaspect([1.5 1 1])
if savefig
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 6];
    fig.PaperPositionMode = 'Manual';
    fig.PaperSize = [6 6];
    print(gcf,  ['../../figures/tuning/population_over_time'],'-dsvg')
end
%%%% plot first figure
colorfig = figure(7); clf;
h = gca;
set(h, 'fontsize',16)
hold on;

n_dv_bins   = numel(results.dv_axis);
n_cells     = size(results.fga_aa_cell,2);
clrs        = color_set(n_dv_bins);
Markersize  = 2;

for i=1:res:n_dv_bins
    plot(h,results.t0s,movmean(results.tuning_population(:,i),3)-0.5,'-','Color',clrs(i,:),'Markerfacecolor',clrs(i,:),'MarkerSize',Markersize);
end
leg_loc                 = 'Northwest';
legend_str{1}           = 'Negative bound';
legend_str{n_dv_bins}   = 'Positive bound';
for i=2:n_dv_bins-1
    legend_str{i} = ['DV = ' num2str(results.dv_axis(i))];
end
xlabel(h,'Time (s)');
ylabel(h,'Mean normalized FR');
title(h,'SVD Population Average');
%ylim([0 1])
ylim([-.5 .5])
cpos = get(colorfig, 'Position');
pbaspect([1.5 1 1])
if savefig
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 6];
    fig.PaperPositionMode = 'Manual';
    fig.PaperSize = [6 6];
    print(gcf,  ['../../figures/tuning/population_over_time_svd'],'-dsvg')
end


%%%% plot second figure
fig = figure(2); clf;
h = gca;
set(h, 'fontsize',16)
hold on;

plot(h,results.dv_axis, results.fga_ta_population(:)-0.5, '-', 'Color', 'k','linewidth',2);
for i=1:res:numel(results.dv_axis)
    hold on;
    eh = errorplot2(h,results.dv_axis(i), results.fga_ta_population(i)-0.5, results.fga_std_ta_pop(i), 'Marker', '', 'Color', clrs(i,:));
end;

delta = mean(diff(results.dv_axis))*2.5;
xlim([results.dv_axis(1)-delta*1.4, results.dv_axis(end)+delta*1.4]);
xticks = [results.dv_axis(1):2:results.dv_axis(end)];
labels = cellfun(@num2str, mat2cell(xticks, ones(size(xticks,1),1),ones(size(xticks,2),1)),'uni',false);
set(h, 'XTick', xticks, 'XTickLabel', labels, 'Box', 'off');
yl = ylim; ylim(yl);
ylim([-0.6 .6])
set(h, 'Ytick', [-0.5 0 0.5], 'yticklabel',{'-0.5','0','0.5'})
xlabel('Accumulation Value (a)');
ylabel('Mean normalized FR');

% plot grand sigmoidal fit
x = results.dv_axis(1):0.1:results.dv_axis(end);
ypred = dyn_sig(results.grand_fit_betas,x);
%plot(h,x,ypred,'m','linewidth',2);
pbaspect([1.5 1 1])
if savefig
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 6];
    fig.PaperPositionMode = 'Manual';
    fig.PaperSize = [6 6];
    print(gcf,  ['../../figures/tuning/population_sigmoid'],'-dsvg')
end




%%% plot third figure
figure(3); clf;
histogram(results.fr_modulation,30,'facecolor','k')
ylabel('Count')
xlabel('Average \Delta FR (Hz)');
set(gca, 'fontsize',16);
pbaspect([1.25 1 1])
if savefig
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 6];
    fig.PaperPositionMode = 'Manual';
    fig.PaperSize = [6 6];
    print(gcf,  ['../../figures/tuning/population_fr_modulation'],'-dsvg')
end

%%% plot fourth figure

threshold = 2;
sigbins = zeros(size(results.t0s));
for i = 1:length(results.cellid)
for j = 1:length(results.t0s)
    frmod = range(results.fga_cell_residual(j,:,i));
    if frmod >= threshold;
        sigbins(j) = sigbins(j) + 1;
    end
end
end
figure(4); clf;hold on
p = sigbins./length(results.cellid);
se = 1.96.*sqrt((p.*(1-p))./length(results.cellid));
plot(results.t0s, 100.*sigbins./length(results.cellid), 'k','linewidth',2)
shadedErrorBar(results.t0s, 100.*p,100.*se,'k')
ylim([ 0 100])
ylabel('% of cells')
title('Cells with \Delta FR > 2 Hz')
xlabel('Time (s)')
set(gca,'fontsize',16);
plot([0.5 0.5], ylim, 'k--')
pbaspect([1.25 1 1])
if savefig
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 6];
    fig.PaperPositionMode = 'Manual';
    fig.PaperSize = [6 6];
    print(gcf,  ['../../figures/tuning/population_fr_modulation_across_time'],'-dsvg')
end



%%% plot fifth figure
threshold = 0.5;
frmods = zeros(length(results.cellid),length(results.t0s));
frmodsn = zeros(length(results.cellid),length(results.t0s));
for i = 1:length(results.cellid)
    for j = 1:length(results.t0s)
        % can also pull from results.frm_time
%        frmod = range(results.fga_cell_residual(j,:,i)); 
        frmod = results.frm_time(j,i);
        frmods(i,j) = frmod;
    end
    frmods(i,:) = movmean(frmods(i,:),3);
    frmodsn(i,:) = movmean(frmods(i,:)./max(frmods(i,:)),3);
end
figure(5); clf; hold on
plot(results.t0s,     frmodsn(results.fr_modulation > threshold,:),'color',[0 0 0]+.5)
plot(results.t0s, mean(frmodsn(results.fr_modulation > 0.5,:)),'k','linewidth',2)
%plot(results.t0s, max(frmodsn(results.fr_modulation > 0.5,:)),'k-')
%plot(results.t0s, min(frmodsn(results.fr_modulation > 0.5,:)),'k')
%shadedErrorBar(results.t0s, mean(frmodsn(results.fr_modulation > 0.5,:)), std(frmodsn(results.fr_modulation >0.5,:)),'k');
xlim([0.25 2])
ylim([0 1])
xlabel('Time (s)')
ylabel('\Delta FR (% of max/cell)')
set(gca,'fontsize',16)
pbaspect([1.25 1 1])
if savefig
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 6];
    fig.PaperPositionMode = 'Manual';
    fig.PaperSize = [6 6];
    print(gcf,  ['../../figures/tuning/fr_modulation_across_time_normalized'],'-dsvg')
end
figure(8); clf; hold on
plot(results.t0s,     results.frm_time(:,results.fr_modulation > threshold),'color',[0 0 0]+.5)
plot(results.t0s, mean(results.frm_time(:,results.fr_modulation > threshold),2),'k','linewidth',2)
%plot(results.t0s, max(frmodsn(results.fr_modulation > 0.5,:)),'k-')
%plot(results.t0s, min(frmodsn(results.fr_modulation > 0.5,:)),'k')
%shadedErrorBar(results.t0s, mean(frmodsn(results.fr_modulation > 0.5,:)), std(frmodsn(results.fr_modulation >0.5,:)),'k');
xlim([0.25 2])
ylim([0 inf])
xlabel('Time (s)')
ylabel('\Delta FR')
set(gca,'fontsize',16)
pbaspect([1.25 1 1])
if savefig
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 6];
    fig.PaperPositionMode = 'Manual';
    fig.PaperSize = [6 6];
    print(gcf,  ['../../figures/tuning/fr_modulation_across_time'],'-dsvg')
end



%%% plot histogram of cell slopes
figure(6); clf; hold on;
histogram(results.nslope_cell(results.cell_dex), 30,'facecolor','k','edgecolor','k');
xlim([0 inf])
plot(nanmedian(results.nslope_cell(results.cell_dex)), max([max(ylim) 15]), 'kv','markerfacecolor','k','markersize',10)
%plot(.145, max([max(ylim) 15]), 'rv','markerfacecolor','r','markersize',10)
%plot(.195, max([max(ylim) 15]), 'bv','markerfacecolor','b','markersize',10)

%ylim([0 20])
set(gca,'fontsize',16)
ylabel('Count')
xlabel('norm slope')
pbaspect([1.25 1 1])
if savefig
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 6];
    fig.PaperPositionMode = 'Manual';
    fig.PaperSize = [6 6];
    print(gcf,  ['../../figures/tuning/population_slope'],'-dsvg')
end

figure(9); clf; hold on;
histogram(results.svd_slope_cell(results.cell_dex), 30,'facecolor','k','edgecolor','k');
xlim([0 inf])
plot(nanmedian(results.svd_slope_cell(results.cell_dex)), max([max(ylim) 15]), 'kv','markerfacecolor','k','markersize',10)
%plot(.145, max([max(ylim) 15]), 'rv','markerfacecolor','r','markersize',10)
%plot(.195, max([max(ylim) 15]), 'bv','markerfacecolor','b','markersize',10)
%ylim([0 20])
set(gca,'fontsize',16)
ylabel('Count')
xlabel('svd slope')
pbaspect([1.25 1 1])
if savefig
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 6];
    fig.PaperPositionMode = 'Manual';
    fig.PaperSize = [6 6];
    print(gcf,  ['../../figures/tuning/population_slope_svd'],'-dsvg')
end

figure(10); clf; hold on;
histogram(results.slope_cell(results.cell_dex), 30,'facecolor','k','edgecolor','k');
xlim([0 inf])
plot(nanmedian(results.slope_cell(results.cell_dex)), max([max(ylim) 15]), 'kv','markerfacecolor','k','markersize',10)
%plot(.145, max([max(ylim) 15]), 'rv','markerfacecolor','r','markersize',10)
%plot(.195, max([max(ylim) 15]), 'bv','markerfacecolor','b','markersize',10)
%ylim([0 20])
set(gca,'fontsize',16)
ylabel('Count')
xlabel('avg slope')
pbaspect([1.25 1 1])
if savefig
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 6];
    fig.PaperPositionMode = 'Manual';
    fig.PaperSize = [6 6];
    print(gcf,  ['../../figures/tuning/population_slope_avg'],'-dsvg')
end





