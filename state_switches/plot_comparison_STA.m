function plot_comparison_STA(savefig)
if nargin < 1
    savefig = 0;
end

% make fig comparing fraction of cells for each alignment
figure(1); clf; hold on;
load(['../figures/STA/fraction_data_model.mat'])
p = movmean(average_siggy,3);
se = 1.96.*sqrt((p.*(1-p))./n);
%shadedErrorBar(time_vec, p.*100, se.*100, 'k',0.1);
plot(time_vec, p.*100, 'k','linewidth',2)

%plot(time_vec, repmat(mean(average_siggy.*100),size(time_vec)),'k--')


load(['../figures/STA/fraction_data_generative.mat'])
p = movmean(average_siggy,3);
se = 1.96.*sqrt((p.*(1-p))./n);
%shadedErrorBar(time_vec, p.*100, se.*100, 'b',0.1);
plot(time_vec, p.*100, 'b','linewidth',2)

%plot(time_vec, repmat(mean(average_siggy.*100),size(time_vec)),'b--')


ylim([0 20])
plot([0 0], [-1 100], 'k--')
plot([.1 .1], [-1 100], 'r--')
set(gca, 'fontsize',16)
ylabel('Significant cells (% of total)')
xlabel(['Time from state switch (s)'])
pbaspect([1 1 1]);
xlim([-.41 .41])
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 6];
fig.PaperPositionMode = 'Manual';
fig.PaperSize = [6 6];
legend({'Model', 'Generative'})
legend boxoff
if savefig
print(fig, ['../figures/STA/fraction_comparison_switch_encoding.svg'],'-dsvg')
end

figure(2); clf; hold on;
load(['../figures/STA/switch_times_model.mat'])
plot(sort(middle_ts),linspace(1,100,length(middle_ts)),'k-','linewidth',2)


load(['../figures/STA/switch_times_generative.mat'])
plot(sort(middle_ts),linspace(1,100,length(middle_ts)),'r-','linewidth',2)
plot([0 0], ylim, 'k--')
ylabel('Cell # sorted by switch time')
xlabel(['Time from state switch (s)'])
set(gca,'fontsize',16)
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 4 4];
fig.PaperPositionMode = 'Manual';
fig.PaperSize = [4 4];
legend({'Model', 'Generative'})
legend boxoff
if savefig
print(fig, ['../figures/STA/distribution_switch_times_comparison.svg'],'-dsvg')
end


