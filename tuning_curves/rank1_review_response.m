dp = set_dyn_path(1);
%%
fht = 2.5;
fw = 1.5.*1.75*fht;
ppos_ = [28 6 fw fht ];
fw = 1.4*fht;
ppos = [0 10 fw fht ];
lw = 2;
aticks      = [-7.5 -2.5 2.5 7.5];
dv_bins = linspace(-6.5,6.5,9);
dvlims  = dv_bins([1 end]);
%%
% load tuning curves for included population
%fig_prefix = 'modelswitch_';
which_switch    = '';
switch_str  = '';
fig_prefix  = '';
fig_type    = '-dsvg';
fig_prefix = [fig_prefix 'zscore_'];
ra_ylab = '\Delta FR (z-score)';
fn = fullfile(dp.model_fits_dir, ['pop_' fig_prefix...
    'tuning_res.mat']);

load(fn,'pop_res')
%% plot VE for ranks 1-5
rank_var = [pop_res(:).rank_var]';
pop_cellids = [pop_res(:).cellid];
[low_ve, low_ve_ind]        = min(rank_var(:,1));
[sorted_ve, ve_sort_ind]    = sort(rank_var(:,1));
nc      = length(sorted_ve);
idx     = [1 round(nc/2) nc];
idx_    = ve_sort_ind(idx);


%%
% plot worst cell


% find(pop_cellids(ve_sort_ind)==18181)
% find(pop_cellids(ve_sort_ind)==16857)
% find(pop_cellids(ve_sort_ind)==18839)
for ii = idx;
    this_ve_ind = ve_sort_ind(ii);
    
    res     = pop_res(this_ve_ind);
    t0s_    = res.t0s;
    alignment   = 'stimstart';
    
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
    
    this_cellid = pop_cellids(this_ve_ind)
    print_fn = @(fh, str, fig_type) print(fh,fullfile(dp.fig_dir,...
        [num2str(this_cellid) '_' str switch_str]),fig_type,'-painters')
    
    
    fh = figure(3); clf
    set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]));
    ax = axes;
    fgta_line_plot(res,'goodtind',~badtind, ...
        'linewidth', lw, 'ax', ax, ...
        'plot_field','fr_given_ta','dvlims',dvlims);
    print_fn(fh,'fgta',fig_type)
    
    
    fh = figure(3); clf
    set(fh,'position',ppos+[5 0 0 0],'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]));
    ax_r = axes;
    fgta_line_plot(res,'goodtind',~badtind,...
        'ax',ax_r,'plot_field','fgta_resid','dvlims',dvlims);
    ylim([-1 1].*max(abs(ylim)))
    
    print_fn(fh,'fgta_resid', fig_type)
    
    fh = figure(4); clf
    set(fh,'position',ppos+[0 -4 0 0],'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]));
    ax = axes;
    fgta_line_plot(res,'goodtind',~badtind,'ax',ax,'plot_field','map_hat','dvlims',dvlims)
    title('rank 1 approximation')
    ylim([-1 1].*max(abs(ylim)))
    text(.25, .9*max(ylim), sprintf('VE=%.1f%%',100*res.rank_var(1)));
    print_fn(fh,'rank1_map',fig_type);

    
    fh = figure(5); clf
    set(fh,'position',ppos+[5 -4 0 0],'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]));
    ax = axes;
    res = compute_rank1_fgta_approx(res);
    rn = 5;
    res.map_hat_5 = (res.u(:,1:rn)'.*res.s(1:rn))'*res.v(:,1:rn)';
    fgta_line_plot(res,'goodtind',~badtind,'ax',ax,...
        'plot_field','map_hat_5','dvlims',dvlims)
    title(sprintf('rank %i approximation',rn))   
    ra_field = 'rank1_ra_n';
    mt_field = 'rank1_mt_n';
    ylim([-1 1].*max(abs(ylim)))
    text(.25, .9*max(ylim), sprintf('VE=%.1f%%',100*res.rank_var(rn)));
    print_fn(fh,'rank5_map',fig_type);
    
    fh = figure(6); clf
    set(fh,'position',ppos + [0 4 0 0],'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]))
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
    print_fn(fh,'rank1_tuning',fig_type);


    mt_ylab = 'FR modulation (z-score)';
    fh = figure(7); clf
    set(fh,'position',ppos + [5 4 0 0],'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]))
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
    if ~isempty(which_switch)
        plot([.5 .5],ylim,':k')
    end
    xlim(x([1 end]))
    
    if ~isempty(which_switch)
        title_str = 'rank 1 $\hat{m}(t-t_c)$';
    else
        title_str = 'rank 1 $\hat{m}(t)$';
    end
    box off
    plot(xlim,[0 0],'k:','linewidth',1)
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
    print_fn(fh,'rank1_mod',fig_type);
    
    
    fh = figure(8); clf
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
    %ylim([-5 5]+[-.5 .5]);
    plot([0 0],ylim,'-k');
    %tfsz = ax_r.Title.FontSize;
    title('$E[\Delta r|a]$','fontsize',14,'fontweight','normal','interpreter','latex')
    %ylim(ax, [-1 1]*7)
    %axt.Position([2 4]) = axpos([2 4]);
    print_fn(fh,'example_fga_tmn',fig_type);
    % plot a psth
    example_cell_psth('cells',res.cellid, 'fpos', ppos + [-4.5 0 0 0])
end
%% test for relationship between # selective timesteps and goodness of rank1
cell_list   = dyn_cells_db;   % That has to be run once to create cell_list
select_str  = 'normmean > 0' ;
cellids     = cell2mat(extracting(cell_list, 'cellid', select_str));
prefp       = cell2mat(extracting(cell_list, 'prefp', select_str));
psr         = cell2mat(extracting(cell_list, 'prefsideright', select_str));
    
%%
cout_auc_file = fullfile(dp.ephys_summary_dir,'cout_auc.mat');
cout_auc = load(cout_auc_file);
good_cells = cout_auc.good_cells;
n_consec_sig = 8;
alpha = .05;
xlim_off = [-1.25 .55];

x0 = xlim_off(1);
        
ind = cout_auc.good_cells;
assert(isequal(cout_auc.cellids(good_cells),pop_cellids(:)))
assert(isequal(cout_auc.cellids(good_cells),cellids(good_cells)))
good_t = cout_auc.cout_t > x0 & cout_auc.cout_t < 0;
sigR  = cout_auc.cout_p(ind,good_t) > 1-alpha/2;
sigL  = cout_auc.cout_p(ind,good_t) < alpha/2;

sigsumL = movsum(sigL, n_consec_sig, 2);
sigsumR = movsum(sigR, n_consec_sig, 2);

sigsumnp  = (sigsumL == n_consec_sig) | (sigsumR == n_consec_sig);

sigsumnp(:,end) = ones(size(sigsumnp,1),1).*n_consec_sig

sigsumnp = sigR .* psr(good_cells) + sigL .* ~psr(good_cells);

nsig = sum(sigsumnp,2);


nsigt = nsig.*mean(diff(cout_auc.cout_t));

fh = figure(10); clf
ax_ = axes;
set(ax_,'position',get(ax,'position'))
this_ppos = ppos+[0 0 -1.25 0]
this_ppos = ppos;
set(fh,'position',this_ppos,'paperposition',[0 0 this_ppos([3 4])],'papersize',ppos_([3 4]));

scatter(nsigt, rank_var(:,1),'o','markerfacecolor',...
    [1 1 1].*.75,'markeredgecolor', [1 1 1].*.75,'linewidth',1.5)
hold on
scatter(nsigt(idx_), rank_var(idx_,1),'o','markerfacecolor',...
    [.65 .25 .65],'markeredgecolor', colors('purple'),'linewidth',1.5)
ylabel('Rank 1 variance explained')
xlabel('Total duration of side-selectivity (s)')

[rho, pval] = corr(nsigt, rank_var(:,1));
text(1, .65, sprintf('\\rho = %.2f', rho),'fontsize',13)

print_fn(fh,'rank1_dur_corr',fig_type);
%%
fh = figure(1); clf
set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],'papersize',ppos_([3 4]));
ax_ = axes;
set(ax_,'position',get(ax,'position'))
histogram(rank_var(:,1),20,'facecolor',[1 1 1].*.5)
hold on
plot([1 1].*mean(rank_var(:,1)),ylim,'k')
ylims = ylim;
plot(rank_var(idx_,1), max(ylims)+ 1, 'v', ...
    'markeredgecolor', colors('purple'), 'markerfacecolor', colors('purple'))
ylim(ylims+[0 1])
%xlim([.45 1.05])

ylabel('# cells')
xlabel('Variance explained by rank 1')
box off
fn = fullfile(dp.fig_dir, 'rank_ve_hist');
print(fh, fn, '-dsvg', '-painters')

fh = figure(2); clf
set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],'papersize',ppos_([3 4]));
ax_ = axes;
set(ax_,'position',get(ax,'position'))
plot(1:size(rank_var,2), rank_var, '-','color',[1 1 1].*.8,'linewidth',1)
hold on
errorbar(1:size(rank_var,2), mean(rank_var), 1.96*sem(rank_var),...
    'ko-','capsize',0,'linewidth',2,'markerfacecolor','w')
box off
ylabel('Variance explained')
xlabel('Rank')
%plot(1:size(rank_var,2), rank_var(idx_,:), 'b')
fn = fullfile(dp.fig_dir, 'rank_ve_lines');
print(fh, fn, '-dsvg', '-painters')
