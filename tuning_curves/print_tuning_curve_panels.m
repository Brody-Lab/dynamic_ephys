krn_width   = 0.1;
force_frdv  = 1;            % keep as one if rebinning
force_bin   = 0;
force_dv    = 1;
norm_type   = 'none';
ops.mask_after_stim_firing    = 1;
ops.mask_stim_firing          = 0;
ops.average_a_vals            = 1;
ops.average_fr                = 1;
ops.min_fr_mod                = 1;
ops.fit_time_sigmoid          = 0;
ops.use_nresidual             = 1;
ops.plot_res                  = 3;
dvlims = [];

ratname = 'H066';
fht = 2.5;
fw = 1.75*fht;
dp          = set_dyn_path;
direction   = 'backward';

krn_type    = 'halfgauss';
end_mask_s  = 0.0;
alignment   = 'stimstart'; 


%% re-analyze cell 18181 
excellid = 18181;
data = dyn_cell_packager(excellid);
[~, vec_data, ~,~] = get_behavior_data(dp.spikes_dir, ...
    data.cellid, data.sessid,ops);
[model, constant_x, allsame] = get_data_model_p(data, vec_data);
[align_strs, align_args] = dyn_align_LUT;
align_ind = strmatch(alignment,align_strs,'exact');
nanfrates   = isnan(data.frate{align_ind});
%%
use_switches = 1
if use_switches
which_switch    = 'model';
switch_params   = struct('t_buffers', [.2 .2]);
else
    which_switch    = [];
    switch_params   = [];
end

do_print        = 1;

dv_bins         = linspace(-6.5,6.5,9);

aticks      = [-7.5 -2.5 2.5 7.5];
lag         = 0.1;            %0.2;
max_dur     = 1;
switch_t0s  = [-.55:.025:.55];
%dv_bins     = [-5.5:2:5.5];
%dv_bins = 6;
switch alignment
    case 'stimend'
        max_dur = 2;
        t0s         = -1-lag:0.025:-0.1;  % 
        t0s         = -max_dur:0.025:-0.1-lag;  % 
    case 'stimstart'
        t0s         = 0.1:0.025:max_dur-lag;  % 
end

min_switch_t = 0;
max_switch_t = 3;
fig_type = '-dsvg';
norm_type = '';
use_fake = 0;
zscore_frates = 0;
demean_frates = 0;
switch_str = '';
do_shuffle_trials = 0;
do_shuffle_fixing_choice = 0;
do_shuffle_switches = 0;
model_ratnoise = [];
model = [];
if ~isempty(which_switch)
    t0s_ = switch_t0s;
    switch_str = ['_' which_switch(1:5) 'switch'];
    badtind = t0s_ > -.20 & t0s_ < .20;
else
    t0s_ = t0s;
    badtind = t0s_ > Inf;
end

if use_fake
    fname = fullfile(dp.model_dir,'simulate_model_ratnoise');
    load(fname, 'asample_ratnoise','model_ratnoise','t');
    fr_gen = -sign(asample_ratnoise);
    %fr_gen = -2*(asample_ratnoise);
    res_fn = @(do_shuffle_switches) dyn_fr_dv_map(18181, 't0s', t0s_, ...
        'frates', fr_gen, 'ft', t,...
        'model', model_ratnoise, 'data', data,...
        'lag', lag, 'alignment', alignment,...
        'krn_width', krn_width, 'krn_type', krn_type,...
        'norm_type', norm_type,'shuffle_trials',0,...
        'n_dv_bins',dv_bins,'end_mask_s',0,...
        'which_switch',which_switch,...
        'switch_params',switch_params,...
        'demean_frates',demean_frates,'zscore_frates',zscore_frates,...
        'do_svd',1, 'shuffle_trials', do_shuffle_trials, ...
        'shuffle_switches', do_shuffle_switches,...
        'min_switch_t',min_switch_t, 'max_switch_t',max_switch_t,...
        'shuffle_trials_fixing_choice',do_shuffle_fixing_choice);

    print_fn = @(fh, str, fig_type) print(fh,fullfile(dp.fig_dir,...
        ['fake_' str switch_str]),fig_type,'-painters');
    
else
    res_fn = @(do_shuffle_switches, model, t0s_) ...
        dyn_fr_dv_map(excellid, 't0s', t0s_, ...
        'data',data, 'model', model, 'lag', lag, ...
        'alignment', alignment,...
        'krn_width', krn_width, 'krn_type', krn_type,...
        'norm_type', norm_type,'shuffle_trials',0,...
        'n_dv_bins',dv_bins,'end_mask_s',end_mask_s,...
        'demean_frates',demean_frates,'zscore_frates',zscore_frates,...
        'which_switch',which_switch, ...
        'switch_params',switch_params,...
        'shuffle_trials', do_shuffle_trials, ...
        'shuffle_switches', do_shuffle_switches,...
        'min_switch_t',min_switch_t, 'max_switch_t',max_switch_t,...
        'shuffle_trials_fixing_choice',do_shuffle_fixing_choice);
    

    print_fn = @(fh, str, fig_type) print(fh,fullfile(dp.fig_dir,[str switch_str]),...
        fig_type,'-painters')
end

res = res_fn(do_shuffle_switches, model, t0s_);
%res.params.switch_params

if ~do_print
    print_fn = @(fh, str, fig_type) fh;
end 

%%
fprintf('\nsvd rank 1 explains %.1f %% of the variance\n',100*res.rank_var(1))
%% Main tuning curve figure
ppos = [1 6 fw fht ];
if ~isempty(which_switch)
    badtind = t0s_ > -.175 & t0s_ < .175;
    xlab = ['time from ' which_switch ' switch (s)'];
    title_str = sprintf('cell %i $E[r|a,t-t_c]$', res.cellid);
else
    badtind = t0s_ > Inf;
    switch alignment
        case 'stimstart'
            xlab  = 'time from stim on (s)';
        case 'stimend'
            xlab = 'time from stim off (s)';
    end
    title_str = sprintf('cell %i $E[r|a,t]$', res.cellid);
end


lw = 2;


dvlims = dv_bins([1 end]);
dvlims = [];

do_normalize_tuning = 1;
if do_normalize_tuning
    ra = res.rank1_ra_n;
    ra = ra - min(ra);
    ra = ra ./ max(ra);
    res.rank1_ra_nn = ra;
    ra_ylab = {'fractional change' 'in firing rate r'};
    mt_ylab = {'modulation (FR)'};
else
    ra_ylab = '\Delta FR';
    mt_ylab = {'modulation' '(fraction of maximum)'};
end
ra_field = 'rank1_ra_nn';
mt_field = 'rank1_mt_n';
ylab = '\Delta FR';
mt = [res.(mt_field)]';
ra = [res.(ra_field)]';

fht = 2.5;
fw = 1.4*fht;
ppos = [0 10 fw fht ];

% example cell E[r(a,t)]
fh = figure(1); clf
set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]));
ax = axes;
[~, cb] = fgta_line_plot(res, 'goodtind',~badtind, ...
    'linewidth', lw, 'ax', ax, ...
    'plot_field','fr_given_ta','dvlims',dvlims);

title(title_str,'fontweight','normal', 'interpreter', 'latex');
xlabel(xlab);
if demean_frates
     ylabel('\Delta FR (spikes/s)');
else
    ylabel('FR (spikes/s)');
    %ylabel('firing rate (spikes/s)')
end
%ylim([39 55])
xlim(res.t0s([1 end]));
axpos = get(ax,'position');
%cb.Location = 'north'
%cb.Position = cb.Position + [.01 .01 -.5 -.02];
title(cb,sprintf('a'));
% cb.Location = 'east'
d = 1/length(res.dv_axis);
cb.Ticks = [0 1] + [1 -1]*(d*1.5);
cb.Position = cb.Position + [0.02 .15 -.02 -.4];
cb.TickLabels = round(10*res.dv_axis([2 end-1]))/10;
%cb.TickLabels = {sprintf('<%.1f', res.dv_axis(2))  sprintf('>%.1f', res.dv_axis(end-1))};

set(ax,'position',axpos);
print_fn(fh,'example_fgta',fig_type);
%% example residual

fh = figure(2); clf
ppos = [0 6 fw fht ];
set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]));
ax_r = axes;
fgta_line_plot(res,'goodtind',~badtind,...
    'ax',ax_r,'plot_field','fgta_resid','dvlims',dvlims);
%title(ax_r,'E[r(a,t)] - E[r(t)]','interpreter','latex')
xlabel(xlab)
ylabel({'\Delta FR (spikes/s)'})
title(ax_r,'$E[\Delta r|a,t]$ for cell 18181','interpreter','latex')

if ~isempty(which_switch)
    title_str = sprintf('cell %i $E[\\Delta r|a,t-t_c]$', res.cellid);
    ylim(ax_r, [-1 1]*10)
else
    title_str = sprintf('cell %i $E[\\Delta r|a,t]$', res.cellid);
    ylim(ax_r, [-1 1]*7)
end
title(title_str,'fontweight','normal','interpreter','latex')

xlim(res.t0s([1 end]))
hold on
plot([0 0],ylim,'-k','linewidth',.1,'color',[1 1 1].*.5)
axpos = get(gca,'position');
print_fn(fh,'example_fgta_resid',fig_type);

%%
fh = figure(4); clf
ppos = [0 2 fw fht ];
set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],...
    'papersize',ppos([3 4]));
axt = axes;
set(axt,'fontsize',12)
fgta_plot_tuning(res,'plot_field','fga_resid_tmn', ...
    'errorbar_field','fga_resid_std','ax',axt,'dvlims',dvlims,...
    'linewidth', 1)
pbaspect(axt, [1 1 1])


xlabel('accumulated value (a)')
ylabel({'\Delta FR (spikes/s)'})
xlim(res.dv_axis([1 end])+[-1 1]);
set(axt,'XTick',aticks);
hold on
ylim([-5 5]+[-.5 .5]);
plot([0 0],ylim,'-k');
tfsz = ax_r.Title.FontSize;
title('$E[\Delta r|a]$','fontsize',tfsz,'fontweight','normal','interpreter','latex')
%ylim(ax, [-1 1]*7)
axt.Position([2 4]) = axpos([2 4]);
print_fn(fh,'example_fga_tmn',fig_type);

%%

res = compute_rank1_fgta_approx(res);

fh = figure(3); clf
ppos = [5 6 fw fht ]
set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]))
ax = axes;
fgta_line_plot(res,'goodtind',~badtind,'ax',ax,'plot_field','map_hat','dvlims',dvlims)
title('rank 1 approximation')
xlabel(xlab)
ylabel({'\Delta FR (spikes/s)'})
axis(ax, 'tight')
ylim(get(ax_r,'ylim'))
print_fn(fh,'example_rank1_map',fig_type);


res.rank1_mt_n2 = res.rank1_mt_n / 2;
res.rank1_ra_n2 = res.rank1_ra_n * 2;
ra_field = 'rank1_ra_n';
mt_field = 'rank1_mt_n';
ra_ylab = {'normalized FR'};
mt_ylab = 'FR modulation (spikes/s)';
%%

fh = figure(5); clf
ppos = [5 2 fw fht ]
set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]))
ax = axes;
fgta_plot_tuning(res,'plot_field',ra_field, ...
    'errorbar_field','','ax',ax,'dvlims',dvlims, 'linewidth', 1)
pbaspect(ax, [1 1 1])
title('rank 1 $\hat{f}(a)$','interpreter','latex')
xlabel('accumulated value (a)')
ylabel(ra_ylab)
box off
ax.TickDir = 'out';
xlim(res.dv_axis([1 end])+[-1 1])
ylim([-.55 .55])
hold on
plot([0 0],ylim,'-k')
set(ax,'XTick',aticks)
print_fn(fh,'example_rank1_tuning',fig_type);
%%
% if ~isempty(which_switch)
%     res_ = res_fn(0, [], t0s);
%     %%
%     fgta_plot_tuning(res_,'plot_field',ra_field, ...
%     'errorbar_field','','ax',ax,'dvlims',dvlims, 'linewidth', 1,...
%     'linecolor',[1 1 1].*.5)
% 
% end
% %%
fh = figure(6); clf
ppos = [10 2 fw fht ]
set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]))
ax = axes;
x = res.t0s;
if ~isempty(badtind)
    xx = res.t0s;
    xx(badtind) = nan;
    yy = res.(mt_field);
    yy(badtind) = nan;
    plot(xx, yy, 'k', 'linewidth', 1.5)
    bad0 = find(badtind,1)-1;
    badn = find(badtind,1,'last')+1;
    hold on
    plot(xx([bad0 badn]), yy([bad0 badn]), ':k', 'linewidth', 1)
    plot([0 0],ylim,':k')
else
    plot(res.t0s, res.(mt_field), 'k', 'linewidth', 1)
    hold on
end
if isempty(which_switch)
    plot([.5 .5],ylim,':k')
end
xlim(x([1 end]))

if ~isempty(which_switch)
title_str = 'rank 1 $\hat{m}(t-t_c)$';
else
    title_str = 'rank 1 $\hat{m}(t)$';
end

    title(title_str,'interpreter','latex')
xlabel(xlab)
ylabel(mt_ylab)
box off
ax.TickDir = 'out';
cb = colorbar
drawnow 
pos = get(ax,'position')
delete(cb)
set(ax,'position',pos)
print_fn(fh,'example_rank1_mod',fig_type);

%% population figure
% load up the "good" cells
cout_auc_file = fullfile(dp.ephys_summary_dir,'cout_auc.mat');
v = load(cout_auc_file)
pop_cellids = v.cellids(v.good_cells);
pop_cellids = pop_cellids(1:end);

%pop_cellids = [18181 17784 ]%16857  ];
pop_sessids = zeros(size(pop_cellids));

for cc = 1:length(pop_cellids)
    pop_sessids(cc) = bdata(['select sessid from cells where ' ...
        'cellid={S}'],pop_cellids(cc));
end
%% population level

refit = 1;
norm_type = 'none';
frbins = [];
demean_frates = 0;
zscore_frates = 1;

switch which_switch
    case ''
        fun = @() pop_dyn_fr_dv_map(pop_cellids, pop_sessids,...
            't0s', t0s, 'lag', lag, ...
            'alignment', alignment,...
            'krn_width', krn_width, 'krn_type', krn_type,...
            'norm_type', norm_type,'shuffle_trials',0,...
            'demean_frates',demean_frates,'zscore_frates',zscore_frates,...
            'n_dv_bins',dv_bins,'end_mask_s',end_mask_s);
        
        fig_prefix = '';
    case 'model'
        fun = @() pop_dyn_fr_dv_map(pop_cellids, pop_sessids,...
            't0s', switch_t0s, 'lag', lag, ...
            'alignment', alignment,...
            'krn_width', krn_width, 'krn_type', krn_type,...
            'norm_type', norm_type,'shuffle_trials',0,...
            'n_dv_bins',dv_bins,'end_mask_s',end_mask_s,...
            'demean_frates',demean_frates,'zscore_frates',zscore_frates,...
            'which_switch','model',...
            'switch_params',switch_params);
        fig_prefix = 'modelswitch_';
    case 'generative'
        fun = @() pop_dyn_fr_dv_map(pop_cellids, pop_sessids,...
            't0s', switch_t0s, 'lag', lag, ...
            'alignment', alignment,...
            'krn_width', krn_width, 'krn_type', krn_type,...
            'norm_type', norm_type,'shuffle_trials',0,...
            'n_dv_bins',dv_bins,'end_mask_s',end_mask_s,...
            'demean_frates',demean_frates,'zscore_frates',zscore_frates,...
            'which_switch','generative',...
            'switch_params',switch_params);
        fig_prefix = 'genswitch_';
end

if zscore_frates
    fig_prefix = [fig_prefix 'zscore_'];
end

if refit | ~exist(fn,'file')
    pop_res = fun();
    keyboard
    if isempty(which_switch)
        fn = fullfile(dp.model_fits_dir, ['pop_' fig_prefix...
            'tuning_res.mat']);
    else
        stadir = get_sta_dirname(pop_res(1).params.switch_params);
        fn = fullfile(stadir, ['pop_' fig_prefix...
            'tuning_res.mat']);
end
    save(fn, 'pop_res','-v7.3')
else
    load(fn,'pop_res')
end

%% plot population level results
%close all
do_normalize_tuning = -0;

fignamefun = @(x) fullfile(dp.fig_dir,[fig_prefix x]);
print_fn = @(fh,x,fig_type) print(fh,fignamefun(x),fig_type,'-painters');


ra_field = 'rank1_ra_n';
mt_field = 'rank1_mt_n';

mt = [pop_res.(mt_field)]';
ra = [pop_res.(ra_field)]';

if do_normalize_tuning
    ra = ra - min(ra,[],2);
    ra = ra ./ max(ra,[],2);
    ylab = {'fractional change' 'in firing rate r'}
else
    ylab = '\Delta FR';
end


dv_axis = pop_res(1).dv_axis;
flip_cell = mean(ra(:,dv_axis>0),2) < mean(ra(:,dv_axis<0),2);

pop_fgta = cat(3,pop_res.fr_given_ta);
pop_fgta_r = cat(3,pop_res.fgta_resid);

ra(flip_cell,:) = fliplr(ra(flip_cell,:));
pop_fgta_r(:,:,flip_cell) = fliplr(pop_fgta_r(:,:,flip_cell)); 
pop_fgta(:,:,flip_cell) = fliplr(pop_fgta(:,:,flip_cell)); 

mean_pop_fgta_r = mean(pop_fgta_r,3);
mean_pop_fgta = mean(pop_fgta,3);
tempres = pop_res(1);
tempres.fgta_resid = mean_pop_fgta_r;
tempres.fr_given_ta = mean_pop_fgta;

tempres = compute_rank1_fgta_approx(tempres,'which_map','fr_given_ta');
tempres.dv_axis = dv_axis;

fh = figure(1); clf
set(fh,'position',[5 5 fw fht])
set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]))

x = pop_res.t0s;
plot(x([1 end]),[0 0],'k')
hold on
plot(x,mt,'color',[1 1 1].*.7)
plot(x,mean(mt),'k','linewidth',2)
%errorbar(x,mean(mt),1.96*nansem(mt),'linewidth',lw,'color','k','capsize',0)
box off

xlabel(xlab)
xlim(x([1 end]))
cb = colorbar
if zscore_frates
ylabel({'modulation' '(z-scores)'})
end

drawnow
pos = get(gca,'position')

delete(cb)
set(gca,'position',pos)
print_fn(fh,['selective_' mt_field],fig_type)

fh = figure(11); clf
set(fh,'position',[5 5 fht fht])
xx = pop_res.dv_axis;
plot(xx,ra,'color',[.5 .5 .5])
hold on
plot(xx,mean(ra),'k','linewidth',2)
%errorbar(xx,mean(ra),1.96*sem(ra),'linewidth',lw,'color','k','capsize',0)

%ylim([-6 6])
box off
ylabel(ylab)
xlabel('accumulated value (a)')
xlim(tempres.dv_axis([1 end])+[-1 1])
if ~do_normalize_tuning & strcmp(ra_field, 'rank1_ra_n')
    ylim([-.61 .61])
elseif do_normalize_tuning
    ylim([-.1 1.1])
end
%ylim([-5.5 5.5])
set(gca, 'TickDir', 'out')
print_fn(fh,['selective_' ra_field],fig_type)
%%
all_var = horzcat(pop_res.rank_var);
all_var = all_var(1,:);
fh = figure(3); clf
set(fh,'position',[5 5 fht fht])
histogram(all_var)
fprintf(['\nfor single cells, svd rank 1 explains an average of '...
    '%.1f %% of the variance SD %.1f%% \n'],100*mean(all_var),100*std(all_var))

%%
fprintf('\n population average svd rank 1 explains %.1f %% of the variance\n',100*tempres.rank_var(1))

%%
pref_clr = [.8 .25 .8];
npref_clr = [.8 .65 .25];
fh = figure(21); clf

set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]))
ax = axes;
fgta_line_plot(tempres,'goodtind',~badtind,'ax',ax,...
    'plot_field','fr_given_ta','dvlims',dvlims,...
    'up_clr', pref_clr, 'down_clr', npref_clr);
ylabel(ax,'\Delta FR (z-score)')
xlabel(ax,xlab)
title('population average','interpreter','latex')
ylim([-1 1].*.175)
print_fn(fh,'pop_fgta',fig_type)
%% population average residual map

fh = figure(20); clf
ppos = [1 1 fw fht ]

set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]))
ax = axes;


[ax, cb] = fgta_line_plot(tempres,'goodtind',~badtind,'ax',ax,...
    'plot_field','fgta_resid','dvlims',dvlims,...
    'up_clr', [.8 .25 .8], 'down_clr', [.8 .65 .25]);
ylabel(ax,'\Delta FR (z-score)')
title('average of selective cells','interpreter','latex')
xlabel(xlab)
if isempty(which_switch)
    ylim([-1 1].*.125)
else
    ylim([-1 1].*.15)
end
title(cb,{'a for ' 'pref. side'})
d = 1/length(tempres.dv_axis)
drawnow
axpos = get(ax,'position')
cb.Ticks = [0 1] + [1 -1]*(d*1.5);
cb.Position = cb.Position + [0.02 .15 -.02 -.4];
cb.TickLabels = round(10*tempres.dv_axis([2 end-1]))/10;

set(ax,'position',get(ax_r,'position'),'ytick',[-.1 0 .1])
print_fn(fh,'pop_fgta_resid',fig_type)


%% POP TUNING + MOD
do_normalize_tuning = 0
this_tmn = tempres.(ra_field);

if do_normalize_tuning
    this_tmn = tempres.(ra_field);
    this_tmn_ = this_tmn - min(this_tmn)
    tempres.(ra_field) = this_tmn_;
end
mt_ylab = 'FR modulation (z-score)';

fh = figure(7); clf
ppos = [5 2 fw fht ]
set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],...
    'papersize',ppos([3 4]))
ax = axes;
fgta_plot_tuning(tempres,'linecolor','k',...
    'plot_field',ra_field, ...
    'errorbar_field','','ax',ax,'dvlims',dvlims,'linewidth',1,...
    'up_clr',pref_clr, 'dwn_clr', npref_clr)
hold(ax,'on')
pbaspect(ax, [1 1 1])
set(ax,'XTick',aticks)
title('rank 1 $\hat{f}(a)$','interpreter','latex')
xlabel('accumulated value (a)')
ylabel(ra_ylab)
ylim([-.55 .55])
plot([0 0],ylim,'-','color',[1 1 1]*.8)
box off
ax.TickDir = 'out';
drawnow
xlim(res.dv_axis([1 end])+[-1 1])
xx = tempres.dv_axis;
this_b = fit_four_param_psycho(xx(:),this_tmn(:));
%plot(xx,fourParamPsychometric(this_b, xx(:))+min(this_tmn(:)),'r');
slope_pop_mean = prod(this_b([2 3]))/4;

print_fn(fh,'pop_ra',fig_type)

%%
if ~isempty(which_switch)
    temp = load(fullfile(dp.model_fits_dir, ['pop_tuning_res.mat']));
    
    %%
    pop_res_ = temp.pop_res;
    pop_fgta_ = cat(3,pop_res_.fr_given_ta);
    pop_fgta_r_ = cat(3,pop_res_.fgta_resid);
    tempres_ = pop_res_(1);
    
    ra_ = [pop_res_.(ra_field)]';
    flip_cell_ = mean(ra_(:,dv_axis>0),2) < mean(ra_(:,dv_axis<0),2);

    
    ra_(flip_cell_,:) = fliplr(ra_(flip_cell_,:));
    
    pop_fgta_r_(:,:,flip_cell_) = fliplr(pop_fgta_r_(:,:,flip_cell_));
    pop_fgta_(:,:,flip_cell_) = fliplr(pop_fgta_(:,:,flip_cell_));
    
    mean_pop_fgta_r_ = mean(pop_fgta_r_,3);
    mean_pop_fgta_ = mean(pop_fgta_,3);
    
    
    tempres_.fgta_resid = mean_pop_fgta_r_;
    tempres_.fr_given_ta = mean_pop_fgta_;
    
    tempres_ = compute_rank1_fgta_approx(tempres_,'which_map','fr_given_ta');
    tempres_.dv_axis = dv_axis;
    fgta_plot_tuning(tempres_,'plot_field',ra_field, ...
    'errorbar_field','','ax',ax,'dvlims',dvlims, 'linewidth', 1,...
    'linecolor',[1 0 1].*.5)
%%
tempflip = res_;
tempflip.(ra_field) = -(tempflip.(ra_field));
fgta_plot_tuning(tempflip,'plot_field',ra_field, ...
    'errorbar_field','','ax',ax,'dvlims',dvlims, 'linewidth', 2,...
    'linecolor',[1 0 1].*.5)
end
%%

fh = figure(8); clf
ppos = [5 5 fw fht ]
set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],...
    'papersize',ppos([3 4]))
ax = axes
hold on
if ~isempty(badtind)
    xx = tempres.t0s;
    xx(badtind) = nan;
    yy = tempres.(mt_field);
    yy(badtind) = nan;
    plot(xx, yy, 'k', 'linewidth', 1.5)
    bad0 = find(badtind,1)-1;
    badn = find(badtind,1,'last')+1;
    hold on
    plot(xx([bad0 badn]), yy([bad0 badn]), ':k', 'linewidth', 1)
    ylim([0 max(ylim)])
    plot([0 0],ylim,'-', 'color', [1 1 1].*.8)
else
    plot(tempres.t0s, tempres.(mt_field), 'k', 'linewidth', 1)
    hold on
    plot([.5 .5],ylim, '-k', 'linewidth', 1)
end
title('rank 1 $\hat{m}(t)$','interpreter','latex')
xlabel(xlab)
box off
ax.TickDir = 'out';
ylim([0 max(ylim)])
xlim(x([1 end]))
ylabel(mt_ylab)
pos = get(ax_r,'position')
set(ax,'position',pos)

print_fn(fh,'pop_frm',fig_type)

%%
fh = figure(11); clf
ppos = [5 5 fw fht ]
set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]))
ax = axes
[~, sortind] = sort(mean(mt,2),'descend');
imagesc(mt(sortind,:),'x',pop_res(1).t0s)
caxis([0 .3])
colormap(bone)
cb = colorbar('northoutside')
title(cb,'modulation (z-scores)')
xlabel('time from model switch')
ylabel('cell #')
print_fn(fh,'pop_frm_heat',fig_type)

%%
mt_n = [pop_res.rank1_mt_n]';
ra_n = [pop_res.rank1_ra_n]';

ra_n(flip_cell,:) = fliplr(ra_n(flip_cell,:));
mean_mt = mean(mt_n,2);
mean_ra_n = mean_mt .* ra_n;
[~, sortorder] = sort(mean_mt,'descend');
figure(111); clf
imagesc(mean_ra_n(sortorder,:),'x',dv_axis)
colormap(colormapRedBlue)
caxis([-2 2]*.1)
colorbar('location','southoutside')
%%
%pop_fga_tmn = vertcat(pop_res.fga_tmn);
fitfun = @(beta, clickdiff) fourParamPsychometric(beta, clickdiff);
p0 = [.1 .8 .01 0]; % [left lapse (right lapse - left lapse) slope controller bias]
opts.MaxIter = 100000;
opts.TolFun = 1e-15;
opts.TolX = 1e-15;
xx = pop_res(1).dv_axis;
slope_field = 'rank1_ra_n';
tmns = horzcat(pop_res.(slope_field))';
figure(1); clf
for ii = 1:length(pop_res)
    %%
    this_tmn = pop_res(ii).(slope_field) ;
    if mean(this_tmn(xx<0)) > mean(this_tmn(xx>0))
        this_tmn = -this_tmn;
    end
    this_tmn = this_tmn - min(this_tmn);
    this_tmn = this_tmn ./ max(this_tmn);
    tmns(ii,:) = this_tmn';

    b(ii,:) = fit_four_param_psycho(xx(:),this_tmn(:))
    
    subplot(121); cla
    plot(xx,fourParamPsychometric(b(ii,:), xx(:)),'r');
    hold on
    plot(xx,this_tmn,'k.-');
       ylim([0 1])
    
    subplot(122);
    
    hold on
    plot(xx,this_tmn,'.-','color',[.5 .5 .5]);
       ylim([0 1])
    drawnow
 
    %%
    %pause()
end
plot(xx,mean(tmns),'k','linewidth',2)

slopes = b(:,3).*b(:,2)/4;
fh = figure(10); clf
ppos = [1 6 fw fht ]
set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]))
histogram(slopes,[0:.025:.325],'facecolor',colors('purple'))
ylabel('number of neurons')
xlabel('slope')
box off
hold on
plot([1 1].*slope_pop_mean,ylim,'k')
%plot([1 1].*mean(slopes),ylim,'color',colors('purple'))
plot(percentile(slopes,[5 95]),[1 1]*max(ylim),'-','color',colors('purple'))
plot(mean(slopes),[1 1]*max(ylim),'.','color',colors('purple'))
print_fn(fh,['pop_slopes' slope_field],fig_type)
%%


%% cross validation

debugcell = 18181;
debugsessid = bdata(['select sessid from cells where cellid={S}'],...
	debugcell); %#ok<NBRAK>

dv_bins = [ -8:.25:8 ];
%dv_bins = 20;

model   = get_data_model_p(debugsessid);

data    = dyn_cell_packager(debugcell);
nt      = length(data.trials.trialnums);
randtrials = randperm(nt);
xval_trials1 = randtrials(1:round(nt/2));
xval_trials2 = randtrials(xval_trials1(end)+1:nt);
res_xval1 = dyn_fr_dv_map(debugcell,...
    't0s', t0s, 'lag', lag, ...
    'frbins', frbins, 'alignment', alignment,...
    'krn_width', krn_width, 'krn_type', krn_type,...
    'norm_type', norm_type,'shuffle_trials',0,...
    'n_dv_bins',dv_bins,'end_mask_s',end_mask_s,...
    'which_trials', xval_trials1,...
    'demean_frates',demean_frates,'zscore_frates',zscore_frates,...
    'data', data);
plot_cell_map(res_xval1, 'fh', 1)

res_xval2 = dyn_fr_dv_map(debugcell,...
    't0s', t0s, 'lag', lag, ...
    'frbins', frbins, 'alignment', alignment,...
    'krn_width', krn_width, 'krn_type', krn_type,...
    'norm_type', norm_type,'shuffle_trials',0,...
    'n_dv_bins',dv_bins,'end_mask_s',end_mask_s,...
    'demean_frates',demean_frates,'zscore_frates',zscore_frates,...
    'which_trials', xval_trials2, 'data', data);
plot_cell_map(res_xval2, 'fh', 2)
%

fh = figure(1); clf
fht = 2.5;
fw = 1.75*fht;
ppos = [1 10 1.5*fw fht ]
dvlims = [-5 5]
lw = 1;
set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]))
ax = subplot(121);
fgta_line_plot(res_xval1,'linewidth',lw,'ax',ax,'plot_field','fr_given_ta','dvlims',dvlims)
title({'E[r(a,t)]' 'Split 1'},'interpreter','latex')
xlabel('time from stim on (s)')
ylabel('FR')
ax = subplot(122);
fgta_line_plot(res_xval2,'linewidth',lw,'ax',ax,'plot_field','fr_given_ta','dvlims',dvlims)
title({'E[r(a,t)]' 'Split 2'},'interpreter','latex')
xlabel('time from stim on (s)')
ylabel('FR')


fh = figure(2); clf
ppos = [1 6 1.5*fw fht ]
set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]))
ax = subplot(121);
fgta_line_plot(res_xval1,'ax',ax,'plot_field',...
    'fgta_resid','dvlims',dvlims)
title({ 'E[r(a,t)] - E[r(t)]' 'Split 1'},'interpreter','latex')
xlabel('time from stim on (s)')
ylabel('\Delta FR')
ylim(ax, [-5 5])
ax = subplot(122);
fgta_line_plot(res_xval2,'ax',ax,'plot_field','fgta_resid','dvlims',dvlims)
title({ 'E[r(a,t)] - E[r(t)]' 'Split 2'},'interpreter','latex')
xlabel('time from stim on (s)')
ylabel('\Delta FR')
ylim(ax, [-5 5])




fh = figure(4); clf
ppos = [1 10 fw fht ]
set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]))
ax = axes;
fgta_plot_tuning(res_xval1,'linecolor','k',...
    'plot_field','fga_tmn', ...
    'errorbar_field','fga_std','ax',ax,'dvlims',dvlims)
hold(ax,'on')
fgta_plot_tuning(res_xval2,'linecolor','r',...
'plot_field','fga_tmn', ...
    'errorbar_field','fga_std','ax',ax,'dvlims',dvlims)
pbaspect(ax, [1 1 1])
title('E[r(a)]','interpreter','latex')
xlabel('accumulated value (a)')
ylabel('\Delta FR')
box off
ax.TickDir = 'out';
%print(fh,fignamefun('example_switch_rank1_tuning'),'-dsvg','-painters')


fh = figure(5); clf
ppos = [5 10 1.5*fw fht ]
set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]))
ax = subplot(131);
fgta_plot_tuning(res_xval1,'linecolor','k',...
    'plot_field','rank1_ra_nm', ...
    'errorbar_field','','ax',ax,'dvlims',dvlims)
hold(ax,'on')
fgta_plot_tuning(res_xval2,'linecolor','r',...
'plot_field','rank1_ra_nm', ...
    'errorbar_field','','ax',ax,'dvlims',dvlims)
%pbaspect(ax, [1 1 1])
title('rank 1 $\hat{f}(a)$','interpreter','latex')
xlabel('accumulated value (a)')
ylabel('\Delta FR')
box off
ax.TickDir = 'out';
ax = subplot(1,3,[2 3]);
plot(res_xval1.t0s, res_xval1.rank1_mt_nm, 'k', 'linewidth', lw)
hold on
plot(res_xval2.t0s, res_xval2.rank1_mt_nm, 'r', 'linewidth', lw)
%pbaspect(ax, [2 1 1])
title('rank 1 $\hat{m}(t)$','interpreter','latex')
xlabel('time from stim on (s)')
ylabel('\Delta FR')
box off
ax.TickDir = 'out';
cb = colorbar
drawnow 
pos = get(ax,'position')
delete(cb)
set(ax,'position',pos)
legend('split 1','split 2','location','southeast')
%print(fh,fignamefun('example_rank1_mod'),'-dsvg','-painters')

%% SWITCH



%%

fh = figure(10); clf

ppos = [5 10 fw fht ]
lw = 2;
set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]))
ax_h = axes;
fgta_heat_plot(res,'linewidth',lw,'ax',ax_h,...
    'plot_field','fgta_resid','dvlims',dvlims)
axis(ax_h, 'tight')
title('$\Delta FR$','interpreter','latex')
xlabel(xlab)
box off
print_fn(fh,'example_fgta_resid_heat',fig_type);

fh = figure(30); clf
ppos = [10 6 fw fht ]
set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]))
ax_map = axes;
fgta_heat_plot(res,'ax',ax_map,'plot_field','map_hat','dvlims',dvlims)
title('rank 1 approximation','interpreter','latex')
xlabel(xlab)
ylabel('\Delta FR')
axis(ax_map, 'tight')
ylim(get(ax_h,'ylim'))
xlim(get(ax_h,'xlim'))
box off
print_fn(fh,'example_rank1_map_heat',fig_type);

fh = figure(7); clf
ppos = [10 10 fw fht ]
set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]))
ax = axes;
plot(1:length(res.rank_var), res.rank_var, 'k', 'linewidth', lw)
pbaspect(ax, [1 2 1])
box off
ax.TickDir = 'out';
ylim([.7 1])
ylabel('variance explained')
xlabel('rank')
xlim([1 length(res.rank_var)])
print_fn(fh,'example_rank_var',fig_type);
%%
switch which_switch
    case 'model'
        badtind = pop_res(1).t0s > -.15 & pop_res(1).t0s < .15;
end
       
%%
ra_field = 'rank1_ra_n';

xlab = '';
fh = figure(1); clf
ppos = [1 6 10 fht ]
set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]))
n_shuffles = 0
ax = subplot(132); hold on
for ss = 1:n_shuffles
    if ~isempty(badtind)
        xx = res.t0s;
        xx(badtind) = nan;
        yy = res.(mt_field);
        yy(badtind) = nan;
        plot(xx, yy, 'color', 'k', 'linewidth', lw)
        bad0 = find(badtind,1)-1;
        badn = find(badtind,1,'last')+1;
        hold on
        plot(xx([bad0 badn]), yy([bad0 badn]), ':k', 'linewidth', lw)
    else
        plot(res.t0s, res.(mt_field), 'k', 'linewidth', lw)
        
    end
end
%pbaspect(ax, [2 1 1])
title('rank 1 $\hat{m}(t)$','interpreter','latex')
xlabel(xlab)
ylabel({'modulation' '(fraction of maximum)'})
box off
ax.YGrid = 'on';
pos2 = get(ax,'position');
linkaxes([ax ax],'x')
plot(xlim,[0 0],'k')
xlim(res.t0s([1 end])+[-.1 .1])
ylim([min(0,min(ylim)) max(ylim)])

ax = subplot(131);
fgta_line_plot(res,'colorbar_location','','ax',ax,...
    'plot_field','map_hat','dvlims',dvlims,'goodtind',~badtind)
title('rank 1 approximation','interpreter','latex')
xlabel(xlab)
ylabel('\Delta FR')
axis(ax, 'tight')
%ylim(get(ax_r,'ylim'))
drawnow
pos1 = get(ax,'position');


ax = subplot(133);
fgta_plot_tuning(res,'plot_field',ra_field, ...
    'errorbar_field','','ax',ax,'dvlims',dvlims, 'linewidth', lw)
pbaspect(ax, [1 1 1])
title('rank 1 $\hat{f}(a)$','interpreter','latex')
xlabel('accumulated value (a)')
ylabel('\Delta FR')
box off
ax.TickDir = 'out';
xlim(res.dv_axis([1 end])+[-1 1])
hold on
plot([0 0],ylim,':k')
ylim([-.55 .55])
%print_fn(fh,'example_rank1_mod',fig_type);
% 
% hold on
% plot(res.dv_axis, ra_shuffle)

%%
figure(10); clf
res.mass_ta = squeeze(sum(res.mass_tfa,2));
res.mass_ta = res.mass_t/10;
res.mass_ta = squeeze(sum(res.mass_tfa,2))./res.mass_t/10;
midp = size(res.mass_ta,2)/2  + [0:1];
subplot(311)
fgta_line_plot(res,'plot_field','mass_ta')
hold on
ylabel('mass')
subplot(312)
fgta_heat_plot(res,'plot_field','mass_ta')
colormap(gca,flipud(bone(200)))
subplot(313)
plot(res.t0s,res.mass_t)

%%

which_switch = 'model';
res_switch = dyn_fr_dv_map(excellid, 't0s', switch_t0s, ...
    'model', model, 'trialnums',1:length(model),...
    'lag', lag, ...
    'frbins', [], 'alignment', alignment,...
    'krn_width', krn_width, 'krn_type', krn_type,...
    'norm_type', norm_type,'shuffle_trials',0,...
    'n_dv_bins',dv_bins,'end_mask_s',0,...
    'which_switch',which_switch,'demean_frates',1);

res_gen = dyn_fr_dv_map(excellid, 't0s', switch_t0s, ...
    'model', model, 'trialnums',1:length(model),...
    'lag', lag, ...
    'frbins', [], 'alignment', alignment,...
    'krn_width', krn_width, 'krn_type', krn_type,...
    'norm_type', norm_type,'shuffle_trials',0,...
    'n_dv_bins',dv_bins,'end_mask_s',0,...
    'which_switch','generative','demean_frates',1);
%%

fh = figure(1); clf
ax = axes;
fgta_line_plot(res_switch,'linewidth',lw,'ax',ax,'plot_field','fr_given_ta','dvlims',dvlims)
title('E[r(a,t)]','interpreter','latex')
ylabel('FR')
xlabel(['time from ' res_switch.which_switch ' switch (s)'])
print(fh,fullfile(dp.fig_dir,'example_modelswitch_fgta'),'-dsvg','-painters')

fh10 = figure(10); clf
set(fh10,'position',get(fh,'position'))
ax = axes;
fgta_line_plot(res_gen,'linewidth',lw,'ax',ax,'plot_field','fr_given_ta','dvlims',dvlims)
title('E[r(a,t)]','interpreter','latex')
ylabel('FR')
xlabel(['time from ' res_gen.which_switch ' switch (s)'])
print(fh10,fullfile(dp.fig_dir,'example_genswitch_fgta'),'-dsvg','-painters')

fh = figure(2); clf
ppos = [1 6 fw fht ]
set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]))
ax = axes;
fgta_line_plot(res_switch,'ax',ax,'plot_field','fgta_resid','dvlims',dvlims)
xlabel(['time from ' res_switch.which_switch ' switch (s)'])
title('E[r(a,t)]-E[r(t)]','interpreter','latex')
ylabel('\Delta FR')
print(fh,fullfile(dp.fig_dir,'example_modelswitch_fgta_resid'),'-dsvg','-painters')

fh20 = figure(20); clf
ppos = [1 6 fw fht ]
set(fh20,'position',ppos,'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]))
ax = axes;
fgta_line_plot(res_gen,'ax',ax,'plot_field','fgta_resid','dvlims',dvlims)
xlabel(['time from ' res_gen.which_switch ' switch (s)'])
title('E[r(a,t)]-E[r(t)]','interpreter','latex')
ylabel('\Delta FR')
print(fh20,fullfile(dp.fig_dir,'example_genswitch_fgta_resid'),'-dsvg','-painters')


fh = figure(3); clf
ppos = [1 6 fw fht ]
set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]))
ax = axes;
fgta_line_plot(res_switch,'ax',ax,'plot_field','map_hat','dvlims',dvlims)
title('rank 1 approximation','interpreter','latex')
xlabel(['time from ' res_switch.which_switch ' switch (s)'])
ylabel('\Delta FR')
print(fh,fullfile(dp.fig_dir,'example_modelswitch_rank1_map'),'-dsvg','-painters')


fh = figure(30); clf
ppos = [1 6 fw fht ]
set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]))
ax = axes;
fgta_line_plot(res_gen,'ax',ax,'plot_field','map_hat','dvlims',dvlims)
title('rank 1 approximation','interpreter','latex')
xlabel(['time from ' res_gen.which_switch ' switch (s)'])
ylabel('\Delta FR')
print(fh,fullfile(dp.fig_dir,'example_genswitch_rank1_map'),'-dsvg','-painters')



fh = figure(4); clf
ppos = [1 2 fw fht ]
set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]))
ax = axes;
fgta_plot_tuning(res_switch,'plot_field','fga_resid_tmn', ...
    'errorbar_field','fga_resid_std','ax',ax,'dvlims',dvlims)
pbaspect(ax, [1 1 1])
title('E[$\Delta$ r(a)]','interpreter','latex')
xlabel('accumulated value (a)')
ylabel('\Delta FR')
print(fh,fullfile(dp.fig_dir,'example_modelswitch_fga_tmn'),'-dsvg','-painters')

fh = figure(40); clf
ppos = [1 2 fw fht ]
set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]))
ax = axes;
fgta_plot_tuning(res_gen,'plot_field','fga_resid_tmn', ...
    'errorbar_field','fga_resid_std','ax',ax,'dvlims',dvlims)
pbaspect(ax, [1 1 1])
title('E[$\Delta$ r(a)]','interpreter','latex')
xlabel('accumulated value (a)')
ylabel('\Delta FR')
print(fh,fullfile(dp.fig_dir,'example_genswitch_fga_tmn'),'-dsvg','-painters')


fh = figure(5); clf
ppos = [5 10 fw fht ]
set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]))
ax = axes;
fgta_plot_tuning(res_switch,'linecolor',dp.model_color,...
    'plot_field','rank1_ra_n', ...
    'errorbar_field','','ax',ax,'dvlims',dvlims)
hold(ax,'on')
fgta_plot_tuning(res_gen,'plot_field','rank1_ra_n', ...
    'errorbar_field','','ax',ax,'dvlims',dvlims)
pbaspect(ax, [1 1 1])
title('rank 1 $\hat{f}(a)$','interpreter','latex')
xlabel('accumulated value (a)')
ylabel('\Delta FR')
box off
ax.TickDir = 'out';
print(fh,fullfile(dp.fig_dir,'example_switch_rank1_tuning'),'-dsvg','-painters')



fh = figure(6); clf
ppos = [5 6 fw fht ]
set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]))
ax = axes;
plot(res_gen.t0s, res_gen.rank1_mt_n, 'k', 'linewidth', lw)
hold on
plot(res_switch.t0s, res_switch.rank1_mt_n, 'color',dp.model_color, 'linewidth', lw)
%pbaspect(ax, [2 1 1])
title('rank 1 $\hat{m}(t)$','interpreter','latex')
%xlabel(['time from ' res_switch.which_switch ' switch (s)'])
xlabel(['time from state switch (s)'])
ylabel('\Delta FR')
box off
ax.TickDir = 'out';
cb = colorbar
drawnow 
pos = get(ax,'position')
delete(cb)
hold on
xlim(ax,[-.1 .1]+res_gen.t0s([1 end]))
plot(xlim, [0 0],'k')
hl = legend('generative','model','location','eastoutside')
hl.Position = hl.Position + [-.05 0 0 0];
title(hl,'switch')
set(ax,'position',pos)

box(hl,'off')
print(fh,fullfile(dp.fig_dir,'example_switch_rank1_mod'),'-dsvg','-painters')

fh = figure(7); clf
ppos = [5 6 fw fht ]
set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]))
ax = axes;
plot(1:length(res_gen.rank_var), res_gen.rank_var, 'k', 'linewidth', lw)
hold on
plot(1:length(res_switch.rank_var), res_switch.rank_var, 'color',dp.model_color, 'linewidth', lw)
pbaspect(ax, [1 2 1])
box off
ax.TickDir = 'out';
ylim([.7 1])
ylabel('variance explained')
xlim([1 length(res_switch.rank_var)])

print(fh,fullfile(dp.fig_dir,'example_switch_rank_var'),'-dsvg','-painters')




