
% plot example trial for explaining analysis for panel A
posterior_example;
%%
lag = .1;
max_t = .75;
min_t = -.75;
% plot example cell STAs for panel B
plotSTA(18181,'lag',lag,'fig_num',1);
%%
plotSTA(16857,'lag',lag,'fig_num',2);
plotSTA(17784,'lag',lag,'fig_num',3);
%% plot population STA for panel C and D
recompute = 0;
cnum = 0;
switch cnum
    case 0
        which_correction_str = '';
    case 1
        which_correction_str = '_bonferroni';
    case 2
        which_correction_str = '_bonferroni_modified';
end
which_correction_str = [which_correction_str  '_0']; 
pop_sta_fn = @(which_switch,lag) plot_population_STA(...
    'which_switch',which_switch,...
    'recompute',recompute,'savefig',1,'correction_num',cnum,'lag',lag,...
    'min_t',min_t,'max_t',max_t);

pop_sta_fn('model',lag)

pop_sta_fn('generative',lag)

%%


[fh ax]= plot_comparison_STA(which_correction_str,1)
ylim(ax,[-.65 65])
print(fh, fullfile(dp.sta_fig_dir,['STA_comparison_' which_correction_str ]),fig_type,'-painters')

%ax.YTick = [0:5:35];
fht = 2.5;
fw = 1.75*fht;
ppos = [8 10 fw fht ]
set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]))
