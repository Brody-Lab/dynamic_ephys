function [fh1 ax] = plot_comparison_STA(which_correction_str,savefig)
fig_type = '-dsvg'

if nargin < 1
    which_correction_str = '_bonferroni_modified_0';
end

if nargin < 2
    savefig = 0;
end

dp = set_dyn_path;
% make fig comparing fraction of cells for each alignment
%%
fh1 = figure(1); clf; hold on;

fht = 2;
fw  = 3.75;
ppos = [fw fht];
% set(fh1,'position',[3 7 ppos],'papersize',[ppos],...
%     'paperposition',[0 0 ppos])
set(fh1,'position',[2 5 fw fht],'papersize', [fw fht])

ns = 1;
load(fullfile(dp.sta_fig_dir,...
    ['fraction_data_generative' which_correction_str '.mat']))
time_gen = time_vec;
p_gen = movmean(average_siggy,ns);
se_gen = sqrt((p_gen.*(1-p_gen))./n);
load(fullfile(dp.sta_fig_dir,...
    ['fraction_data_model' which_correction_str '.mat']))
p_mod = movmean(average_siggy,ns);
se_mod = sqrt((p_mod.*(1-p_mod))./n);
time_mod = time_vec;

gen_color = [.5 .5 .5];

patch([-1 -1 1 1].*.55, [0 1 1 0]*5,[1 1 1].*.95, 'edgecolor',[1 1 1].*.95)
hold on
plot([0 0], [-1 100], 'k:','linewidth',1)
shadedErrorBar(time_gen, p_gen.*100, se_gen.*100,{'color',gen_color},0.01);

shadedErrorBar(time_mod, p_mod.*100, se_mod.*100, {'color',dp.model_color},0.1);

ylabel('significant cells (% of total)')
xlabel(['time from state switch (s)'])
ylim([-.1 50])
drawnow
%hl = legend({'Generative', 'Model'},'fontsize',dp.fsz,'location','west')
text(-.45, .98*max(ylim), 'triggered on','fontsize',dp.fsz)
text(-.45, .98*max(ylim)-.1*diff(ylim), 'generative switches','color',gen_color,'fontsize',dp.fsz)
text(-.45, .98*max(ylim)-.2*diff(ylim), 'model switches','color',dp.model_color,'fontsize',dp.fsz)
ax = gca;
xlim(ax,[-.55 .55])
%plot([0 0], [-1 100], 'k--')
%plot(xlim,5*[1 1],'k:')
%%
legend boxoff
if savefig
    print(fh1, fullfile(dp.fig_dir,...
        ['STA_comparison_' which_correction_str ]),fig_type,'-painters')
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
print(fh, fullfile(dp.sta_fig_dir, 'distribution_switch_times_comparison.svg'),fig_type)
end


