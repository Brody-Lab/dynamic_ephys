% script for running all cells through the switch triggered average shuffle
% test, or just loading the files if they already exist

which_switch = 'generative';
[res, dprime, pvals,cellids] = load_STA_population(which_switch, 0,1);

figure(1); clf
pval_plot_lags = res{1}.lags > -.5 & res{1}.lags < 1;
clrs = color_set(6);
posdex = pval_plot_lags > -0.0001;
negdex = pval_plot_lags <  0.0001;

%cmap = dyn_cmap(20).^.6;
cmap = [repmat(clrs(1,:),10,1); repmat([1 1 1], 200, 1); repmat(clrs(end,:),10,1)];
imagesc(sort_by_peak(pvals),'x',res{1}.lags(pval_plot_lags),[0 1]); 
temp = abs(pvals - 0.5);
temp(temp < 0.45) = 0;
sort_dex = sum(temp,2);
sort_dex(isnan(sort_dex)) = 0;
[~, sort_ind] = sort(sort_dex,1 ,'descend');
imagesc(pvals(sort_ind,:),'x',res{1}.lags(pval_plot_lags),[0 1]); 
hold on;
plot([0 0],ylim,'k','linewidth',2)
colormap(cmap)
cb = colorbar;
set(gca,'fontsize',18)
set(gcf,'color','w')
ylabel('Neuron # (sorted by encoding strength)')
xlabel('Time from model switch (s)')
%title({'state switch encoding'});% ['for neurons with ' select_str]})
ylabel(cb, 'd'' percentile w/in shuffled trials')
pbaspect([1.75 2 1])
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8 8];
fig.PaperPositionMode = 'Manual';
fig.PaperSize = [8 8];
print(fig, ['../figures/STA/' which_switch 'switch_encoding.svg'],'-dsvg')

% look at fraction of significant cells
bin_temp = temp;
bin_temp(temp > 0) = 1;
average_siggy = sum(bin_temp,1)./size(temp,1);
figure(2); clf;
time_vec = res{1}.lags(pval_plot_lags);
plot(time_vec, average_siggy.*100, 'k','linewidth',2)
ylim([0 25])
hold on
plot([0 0], [-1 100], 'k--')
set(gca, 'fontsize',16)
ylabel('Significant cells (% of total)')
xlabel(['Time from ' which_switch ' switch (s)'])
pbaspect([1 2 1]);
xlim([-.41 .41])
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 6];
fig.PaperPositionMode = 'Manual';
fig.PaperSize = [6 6];
print(fig, ['../figures/STA/fraction_' which_switch 'switch_encoding.svg'],'-dsvg')

% Save data from this analysis
save(['../figures/STA/fraction_data_' which_switch '.mat'], 'time_vec','average_siggy','pvals')

if 0
% make fig comparing fraction of cells for each alignment
figure(3); clf;
load(['../figures/STA/fraction_data_model.mat'])
plot(time_vec, average_siggy.*100, 'k','linewidth',2)
hold on;
load(['../figures/STA/fraction_data_generative.mat'])
plot(time_vec, average_siggy.*100, 'r','linewidth',2)
ylim([0 25])
plot([0 0], [-1 100], 'k--')
set(gca, 'fontsize',16)
ylabel('Significant cells (% of total)')
xlabel(['Time from state switch (s)'])
pbaspect([1 2 1]);
xlim([-.41 .41])
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 6];
fig.PaperPositionMode = 'Manual';
fig.PaperSize = [6 6];
legend({'Model', 'Generative'})
legend boxoff
print(fig, ['../figures/STA/fraction_comparison_switch_encoding.svg'],'-dsvg')
end


