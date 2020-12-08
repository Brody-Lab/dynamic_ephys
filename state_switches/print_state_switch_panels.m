
% plot example trial for explaining analysis for panel A
posterior_example;
% plot example cell STAs for panel B
plotSTA(18181);
plotSTA(16857);
plotSTA(17784);
%% plot population STA for panel C and D
plot_population_STA('model',1,1);
plot_population_STA('generative',1,1);
%%
which_correction_str = '_bonferroni_modified_0';

plot_comparison_STA(which_correction_str,1)