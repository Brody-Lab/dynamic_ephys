check_switch_inclusion;
%%
strength_histograms;
%% plot example trial for explaining analysis for panel A
posterior_example;
%%
plotSTA_example_switches; 
%%
lag = .1;
max_t = .55;
min_t = -.55;
which_switch = 'model';
t_buffers = [.2 .2];
plot_sta_fn = @(cellid,fig_num) plotSTA(cellid,'which_switch',which_switch,...
    'lag',lag,'fig_num',fig_num,'recompute',0,'ylims',[-1 1].*11,...
    't_buffers', t_buffers, 'alpha',.05)
% plot example cell STAs for panel B
[res, fh, axsta] = plot_sta_fn(18181,1);

title('example cell STR')
fht = 2;
fw  = 3.75;
set(fh,'position',[2 5 fw fht],'papersize', [fw fht])
fig_name = [ num2str(res.cellid) '_' which_switch '_STA'];
print(gcf,  fullfile(dp.fig_dir,fig_name),'-dsvg','-painters')
%%
% plot_sta_fn(16857,2);
% plot_sta_fn(17784,3);
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



pop_sta_fn = @(which_switch,lag,recompute) plot_population_STA(...
    'which_switch',which_switch,...
    'recompute',recompute,'savefig',1,'correction_num',cnum,'lag',lag,...
    'min_t',min_t,'max_t',max_t,...
    'bad_strength', 0,...
    't_buffers', t_buffers);
%%

axpopsta = pop_sta_fn('model',lag,recompute)

%%
recompute = 0;
axpopsta = pop_sta_fn('generative',lag,recompute)

%%
[fh, ax] = plot_comparison_STA(which_correction_str,1);
xlim(ax,[-.55 .55])
ylim(ax,[-.65 50])
ylabel(ax,'% significant')
set(fh,'position',get(figure(1),'position'),...
    'papersize',get(figure(1),'papersize'),...
    'paperposition',get(figure(1),'paperposition'))
%set(ax,'position',get(axpopsta,'position'))
%set(ax,'position',[2 5 fw fht])
print(fh, fullfile(dp.fig_dir,['STA_comparison_' which_correction_str ]),fig_type,'-painters')
