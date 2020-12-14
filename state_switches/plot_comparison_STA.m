function [fh1 ax] = plot_comparison_STA(which_correction_str,savefig)
if nargin < 1
    which_correction_str = '_bonferroni_modified_0';
end

if nargin < 2
    savefig = 0;
end

dp = set_dyn_path;
% make fig comparing fraction of cells for each alignment

fh1 = figure(1); clf; hold on;
ns = 1;
load(fullfile(dp.sta_fig_dir,['fraction_data_generative' which_correction_str '.mat']))
p_gen = movmean(average_siggy,ns);
se_gen = 1.96.*sqrt((p_gen.*(1-p_gen))./n);
load(fullfile(dp.sta_fig_dir, ['fraction_data_model' which_correction_str '.mat']))
p_mod = movmean(average_siggy,ns);
se_mod = 1.96.*sqrt((p_mod.*(1-p_mod))./n);

gen_color = [.5 .5 .5];
%plot(time_vec, p_gen.*100, 'k','linewidth',2)
%plot(time_vec, p_mod.*100, 'color',dp.model_color,'linewidth',2)
%plot(time_vec, repmat(mean(average_siggy.*100),size(time_vec)),'b--')
shadedErrorBar(time_vec, p_gen.*100, se_gen.*100,{'color',gen_color},0.01);

shadedErrorBar(time_vec, p_mod.*100, se_mod.*100, {'color',dp.model_color},0.1);

%plot(time_vec, repmat(mean(average_siggy.*100),size(time_vec)),'k--')
ppos = [8 10 6 3];
set(fh1,'position',ppos,'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]))

%ylim([0 10])
plot([0 0], [-1 100], 'k--')
%plot([.1 .1], [-1 100], 'r--')
ylabel('significant cells (% of total)')
xlabel(['time from state switch (s)'])
%pbaspect([1 1 1]);
xlim([-.41 .41])

%hl = legend({'Generative', 'Model'},'fontsize',dp.fsz,'location','west')
text(-.4, 29, 'triggered on','fontsize',dp.fsz)
text(-.4, 26, 'generative switches','color',gen_color,'fontsize',dp.fsz)
text(-.4, 23, 'model switches','color',dp.model_color,'fontsize',dp.fsz)
ax = gca;
ax.TickDir = 'out'
xlim(ax,[-.45 .75])

legend boxoff
if savefig
print(fh1, fullfile(dp.sta_fig_dir,['STA_comparison_' which_correction_str '.svg']),'-dsvg','-painters')
end

fh = figure(2); clf; hold on;


load(fullfile(dp.sta_fig_dir, ['switch_times_generative' which_correction_str '.mat']))
plot(sort(middle_ts),linspace(1,100,length(middle_ts)),'k-','linewidth',1)
load(fullfile(dp.sta_fig_dir, ['switch_times_model' which_correction_str '.mat']))
plot(sort(middle_ts),linspace(1,100,length(middle_ts)),'-','color',dp.model_color,'linewidth',1)
plot([0 0], ylim, 'k--')
ylabel('Cell # (sorted by switch time)')
xlabel(['time from state switch (s)'])
ppos = [12 10 6 3 ];

set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]))

legend({ 'Generative', 'Model'})
legend boxoff
if savefig
print(fh, fullfile(dp.sta_fig_dir, 'distribution_switch_times_comparison.svg'),'-dsvg')
end


