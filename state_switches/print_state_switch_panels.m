
% plot example trial for explaining analysis for panel A
posterior_example;
% plot example cell STAs for panel B
plotSTA(18181);
plotSTA(16857);
plotSTA(17784);
%% plot population STA for panel C and D
cnum = 2;
switch cnum
    case 0
        which_correction_str = '';
    case 1
        which_correction_str = '_bonferroni';
    case 2
        which_correction_str = '_bonferroni_modified';
end
which_correction_str = [which_correction_str  '_0']; 
plot_population_STA('model',1,cnum);
plot_population_STA('generative',1,cnum);
%%


[fh ax]= plot_comparison_STA(which_correction_str,1)
ylim(ax,[-.5 35])
print(fh, fullfile(dp.sta_fig_dir,['STA_comparison_' which_correction_str '.svg']),'-dsvg','-painters')

%ax.YTick = [0:5:35];
fht = 2.5;
fw = 1.75*fht;
ppos = [8 10 fw fht ]
set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]))
