function [fh1, ax, res] = plot_comparison_STA(which_correction_str,savefig)
fig_type = '-dsvg';

if nargin < 1
    which_correction_str = '_bonferroni_modified_0';
end

if nargin < 2
    savefig = 0;
end

dp = set_dyn_path;
% make fig comparing fraction of cells for each alignment
%%
fh1     = figure(11); clf; hold on;
fht     = 2;
fw      = 3.75;
ppos    = [fw fht];
ns      = 1;
gen_clr = [.5 .5 .5];
set(fh1,'position',[2 5 fw fht],'papersize', [fw fht]);

load(fullfile(dp.sta_dir,...
    ['fraction_data_generative' which_correction_str '.mat']))
time_gen = time_vec;
p_gen = movmean(average_siggy,ns);
se_gen = sqrt((p_gen.*(1-p_gen))./n);

load(fullfile(dp.sta_dir,...
    ['fraction_data_model' which_correction_str '.mat']))
p_mod = movmean(average_siggy,ns);
se_mod = sqrt((p_mod.*(1-p_mod))./n);
time_mod = time_vec;

patch([-1 -1 1 1].*.55, [0 1 1 0]*5,[1 1 1].*.95, 'edgecolor',[1 1 1].*.95)
hold on
plot([0 0], [-1 100], 'k:','linewidth',1)
shadedErrorBar(time_gen, p_gen.*100, se_gen.*100,{'color',gen_clr},0.01);

shadedErrorBar(time_mod, p_mod.*100, se_mod.*100, {'color',dp.model_color},0.1);

ylabel('significant cells (% of total)')
xlabel(['time from state switch (s)'])
ylim([-.1 50])
drawnow
text(-.45, .98*max(ylim), 'triggered on','fontsize',dp.fsz)
text(-.45, .98*max(ylim)-.1*diff(ylim), 'generative switches','color',gen_clr,'fontsize',dp.fsz)
text(-.45, .98*max(ylim)-.2*diff(ylim), 'model switches','color',dp.model_color,'fontsize',dp.fsz)
ax = gca;
xlim(ax,[-.55 .55])

res.model_summary = struct('frac_significant_mn', p_mod*100, ...
    'frac_significant_sem', se_mod*100, 'time', time_mod);
res.gen_summary = struct('frac_significant_mn', p_gen*100, ...
    'frac_significant_sem', se_gen*100, 'time', time_gen);


