function plot_population_STA(which_switch, savefig, correction, bad_strength) 

if nargin < 1
    which_switch = 'model';
end
if nargin < 2
    savefig = 1;
end
if nargin < 3
    correction = 2;
end
if nargin < 4
    bad_strength = 0;
end

force = 0;

dp = set_dyn_path;

[res, dprime, ptiles,cellids, computed,not_computed] = ...
    load_STA_population(which_switch, force,1, bad_strength);

%%
ncells      = length(cellids);
nlags       = length(res{1}.lags);
leftSTA     = nan(ncells,nlags);
rightSTA    = nan(ncells,nlags);
for cc = 1:length(cellids)
    leftSTA(cc,:)   = res{cc}.STR_left_real;
    rightSTA(cc,:)  = res{cc}.STR_right_real;
end
%%
cout_auc_file = fullfile(dp.ephys_summary_dir,'cout_auc.mat');
cf = load(cout_auc_file,'cellids','good_cells');
good_ind = ismember(cellids,cf.cellids(cf.good_cells));




diffSTA   = leftSTA - rightSTA;

which_cells = 1:ncells;
%which_cells = [21  510 474];
fh = figure(10); clf
set(fh, 'position',[5 5 3 6],'papersize',[3 3])
ax = subplot(211);
diffSTA = diffSTA(which_cells,:);
NA = ~isnan(diffSTA);
imagesc(diffSTA,'x',res{1}.lags,'Alphadata',NA,'parent',ax)
set(ax,'color',[1 1 1].*.8)

caxis([-2 2])
cm = flipud(colormapLinear(dp.left_color,50));
cm = [cm; colormapLinear(dp.right_color,50)];
colormap(colormapRedBlue)
colormap(cm)
xlim([-.5 .75])
hold on
plot([0 0],ylim,'k--')
%%
ax2 = subplot(212);
errorbar(res{1}.lags, nanmean(abs(diffSTA(:,:))),nansem(abs(diffSTA(:,:))));
xlim([-.5 .75])
%%

% make map of time points for which each cell is encoding signficantly
%figure(1); clf
figure();
pval_plot_lags = res{1}.lags > -.5 & res{1}.lags < 1;
clrs = color_set(6);
posdex = (res{1}.lags > -0.0001);
negdex = (res{1}.lags <  0.0001);
posdex = posdex(pval_plot_lags);
negdex = negdex(pval_plot_lags);
cmap = [repmat(clrs(1,:),10,1); repmat([1 1 1], 200, 1); repmat(clrs(end,:),10,1)];

if correction == 1 % bonferroni corrections
    m = numel(ptiles);
    threshold = 0.5 - 0.05/ m;
    correction_str = '_bonferroni';
elseif correction == 2 % modified bonferroni, because tests are correlated
    % threshold determined by emperically checking for false positives in the shuffle data, and asking what index keeps it below 0.05
    m = 5;
    threshold = 0.5 - 0.05/ m;
    correction_str = '_bonferroni_modified';
else
    threshold = 0.5 - 0.05;
    correction_str = '';
end

correction_str = [correction_str '_' num2str(bad_strength)];
plot_map = ptiles;
plot_map(:,posdex) = 1-plot_map(:,posdex);
if 1
    plot_map(isnan(plot_map)) = 0.5;
    temp = plot_map - 0.5;
    temp(temp < threshold & temp > -threshold) = 0;
    plot_map = temp + 0.5;
else
    pvals(isnan(pvals)) = 0.5;
    temp = abs(pvals - 0.5);
    temp(temp < threshold) = 0;
end
sort_dex = sum(temp,2);
[~, sort_ind] = sort(sort_dex,1,'descend');

plot_sig_only=false;
if plot_sig_only
    sorted_plot_map = plot_map(sort_ind,:);
    all_zero_dex = all(sorted_plot_map == 0.5,2);
    sorted_plot_map(all_zero_dex,:) = [];
    imagesc(sorted_plot_map,'x',res{1}.lags(pval_plot_lags),[0 1]); 
else
    imagesc(plot_map(sort_ind,:),'x',res{1}.lags(pval_plot_lags),[0 1]); 
end

hold on;
plot([0 0],ylim,'k','linewidth',2)
colormap(cmap)
cb = colorbar;
set(gca,'fontsize',18)
set(gcf,'color','w')
ylabel('Neuron # (sorted by encoding strength)')
xlabel(['Time from ' which_switch ' state change (s)'])
ylabel(cb, 'd'' percentile w/in shuffled trials')
pbaspect([1.75 2 1])
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8 8];
fig.PaperPositionMode = 'Manual';
fig.PaperSize = [8 8];

if savefig
    set(gcf,'renderer','opengl')
    fn = fullfile(dp.sta_fig_dir, [which_switch 'switch_encoding' correction_str]);
    print(fig, fn,'-depsc','-opengl')
end

%%%%%%%
if strcmp(which_switch,'model')
    % this sorted cells file is created by orig_dyn_sequence_psth
    load(fullfile(dp.ephys_summary_dir, 'sorted_cells.mat'));
    % sorted_cellids, sorted_post_stim_cells, sorted_stim_cells
    cellids_sorted = cellids(sort_ind);
    for i=1:length(sorted_post_stim_cells)
        dex = find(sorted_post_stim_cells(i) == cellids_sorted);
        if length(dex) == 1
        plot(0.05,dex,'mo','markerfacecolor','m','markersize',1)
        end
    end
    for i=1:length(sorted_stim_cells)
        dex = find(sorted_stim_cells(i) == cellids_sorted);
        plot(0.1,dex,'bo','markerfacecolor','b','markersize',1)
    end
    if savefig
        set(gcf,'renderer','opengl')
        fn = fullfile(dp.sta_fig_dir, [which_switch 'switch_encoding' correction_str '_sequence_labels']);
        print(fig, fn,'-depsc','-opengl')
    end
end


% look at fraction of significant cells
temp                = abs(temp(good_ind,:));
bin_temp            = temp;
bin_temp(temp > 0)  = 1;
average_siggy       = sum(bin_temp,1)./size(temp,1);
n = size(temp,1);
p = movmean(average_siggy,3);
se = 1.96.*sqrt((p.*(1-p))./n);
time_vec            = res{1}.lags(pval_plot_lags);
%%
figure(2); clf;
ax = axes;
shadedErrorBar(time_vec, average_siggy.*100, se.*100, 'k')
%plot(time_vec, movmean(average_siggy.*100,3), 'k','linewidth',2)
ylim([0 20])
hold on
plot([.1 .1], ylim, 'r--')
plot([0 0], [-1 100], 'k--')
ylabel('Significant cells (% of total)')
xlabel(['time from ' which_switch ' switch (s)'])
pbaspect([1 1 1]);
xlim([-.45 .75])
ax.TickDir = 'out'
box(ax,'off')
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 3 3];
fig.PaperPositionMode = 'Manual';
fig.PaperSize = [3 3];

if savefig
    figsavefname = ['fraction_' which_switch 'switch_encoding' correction_str '.svg'];
    print(fig, fullfile(dp.fig_dir, figsavefname) ,'-dsvg','-painters')
    datasavefname = ['fraction_data_' which_switch correction_str '.mat'];
    % Save data from this analysis
    save(fullfile(dp.sta_fig_dir, datasavefname), 'time_vec','average_siggy','n')
end
%%
%%%% look at distribution of dprime values
figure(3); clf;
ax = axes;
dabs = abs(dprime);
mean_dabs = nanmean(dabs,2);
histogram(mean_dabs(mean_dabs < 10),50, 'facecolor','k')
ylim([0 inf])%
xlim([0 .5])
ylabel('count')
xlabel('mean d'' per cell')
ax.TickDir = 'out';
box(ax,'off')
pbaspect([1 1 1])
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 3 3];
fig.PaperPositionMode = 'Manual';
fig.PaperSize = [3 3];
if savefig
    figsavefname = ['distribution_dprime_' which_switch '.svg'];
    print(fig, fullfile(dp.sta_fig_dir, figsavefname),'-dsvg')
end
%%
%%%% look for consistent encoding cells
const_map   = ptiles;
const_map(const_map > 0.99)     = 1;
const_map(const_map < 0.01)     = -1;
const_map(const_map < 1 & const_map > -1) = 0;
const_pmap   = plot_map;
const_pmap(const_pmap > 0.95)     = 1;
const_pmap(const_pmap < 0.05)     = -1;
const_pmap(const_pmap < 1 & const_pmap > -1) = 0;
const_dex   = zeros(size(const_map));
Q1          = zeros(size(const_map));
Q2          = zeros(size(const_map));
Q3          = zeros(size(const_map));
Q4          = zeros(size(const_map));
isconst     = zeros(1,size(const_map,1));
tconst      = NaN(1,size(const_map,1));
tlastconst  = NaN(1,size(const_map,1));
for i = 1:size(const_map,1)
    first   = -inf;
    last    = inf;
    thisc   = false;

    has_sig         = sum(const_pmap(i,:)~=0) >0;
    const_encode    = sum(const_pmap(i,:)~=0) == 1 || ~any(diff(const_pmap(i,const_pmap(i,:)~=0)));
    one_map_flip    = sum(diff(const_map(i,const_map(i,:)~=0))~=0)==1;
    no_map_flip     = sum(diff(const_map(i,const_map(i,:)~=0))~=0)==0;

    if ~has_sig;
        % this cell never encodes anything
        thisc = false;
        first = -inf;
        last = inf;
    elseif const_encode
    % this cell is always consistent, clean up first and last
    % first should be last significant index before t=0, or t=-.5 if none before
    % last should be first significant index after t=0, or t = 1 if none after
        last =  find(const_map(i,:)~=0 & posdex,1,'first')-1;
        first = find(const_map(i,:)~=0 & negdex,1,'last')+1;
        if isempty(first); first = 1; end
        if isempty(last); last = 59;end
        thisc = true;
    elseif one_map_flip %% What does this case mean???
        % this cell changes its encoding once
        thisc = true;
        switchdex = find(diff(const_map(i,const_map(i,:)~=0)));
        possibledex = find(const_map(i,:)~=0,switchdex+1);
        first = possibledex(switchdex)+1; 
        last  = possibledex(switchdex+1)-1;
        if first-last == 1; last = first; end;
        if first > last; error('bad index'); end;

    elseif no_map_flip %% What does this case mean???
        thisc = true;
        first = find(const_map(i,:)~=0 & posdex,1,'last')+1;
        last = 59;
        if isempty(first); first = find(posdex,1); end
        if first > length(posdex); first = last; end
    else   
        % this cell changes its encoding
    end; 

    isconst(i)      = thisc;
    tconst(i)       = first;
    tlastconst(i)   = last;
end

% get indexes of consistent encoding cells
isconst         = logical(isconst);
cellnums        = 1:length(tconst);
isconst_sort    = isconst(sort_ind);
cs_sort         = cellnums(logical(isconst_sort));
notcs_sort      = cellnums(~logical(isconst_sort));
cs              = cellnums(logical(isconst));
notcs           = cellnums(~logical(isconst));

tconst(isinf(tconst))           = NaN;
tconst_sort     = tconst(sort_ind);
ts              = time_vec(tconst(~isnan(tconst)));
ts_sort         = time_vec(tconst_sort(~isnan(tconst_sort)));

tlastconst(isinf(tlastconst))   = NaN;
tlastconst_sort = tlastconst(sort_ind);
tslast          = time_vec(tlastconst(~isnan(tlastconst)));
tslast_sort     = time_vec(tlastconst_sort(~isnan(tlastconst_sort)));

middle_ts       = mean([ts; tslast]);
middle_ts_sort       = mean([ts_sort; tslast_sort]);

figure(4); clf; hold on;
sortplot = true;
plot_map(plot_map > 0.95) =1;
plot_map(plot_map < 0.05) =0;
if sortplot
    imagesc(plot_map(sort_ind,:),'x',res{1}.lags(pval_plot_lags),[0 1]); 
%    sorted_plot_map = plot_map(sort_ind,:);
%    all_zero_dex = all(sorted_plot_map == 0.5,2);
%    sorted_plot_map(all_zero_dex,:) = [];
%    imagesc(sorted_plot_map,'x',res{1}.lags(pval_plot_lags),[0 1]); 
else
    imagesc(plot_map,'x', res{1}.lags(pval_plot_lags),[0 1]); colormap(cmap);
end
hold on;
plot([0 0],ylim,'k','linewidth',2)
colormap(cmap)
cb = colorbar;
set(gca,'fontsize',18)
set(gcf,'color','w')
ylabel('Neuron # (sorted by encoding strength)')
xlabel(['Time from ' which_switch ' state change (s)'])
ylabel(cb, 'd'' percentile w/in shuffled trials')
pbaspect([1.75 2 1])
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8 8];
fig.PaperPositionMode = 'Manual';
fig.PaperSize = [8 8];
if sortplot
%    plot(ts_sort, cs_sort,'ro','markerfacecolor','r');
%    plot(tslast_sort, cs_sort,'ro')
    plot(middle_ts_sort,cs_sort,'mo','markerfacecolor','m')
%    plot([-0.5 1], [notcs_sort; notcs_sort], 'k-')
else
%    plot(ts, cs,'ro','markerfacecolor','r')
%    plot(tslast, cs,'ro')
    plot(middle_ts,cs,'mo','markerfacecolor','m')
%    plot([-0.5 1], [notcs; notcs], 'k-')
end
numrows = size(plot_map,1);
ylim([0 numrows+1])
if savefig
    fn = [which_switch '_switch_labeled_switch_encoding' correction_str];
    print(fig, fullfile(dp.sta_fig_dir, fn),'-depsc','-opengl')
end


figure(5); clf; hold on;
plot(sort(middle_ts),1:length(middle_ts),'k-','linewidth',2)
plot([0 0], ylim, 'k--')
ylabel('Cell # sorted by switch time')
xlabel(['Time from ' which_switch ' switch (s)'])
set(gca,'fontsize',16)
pbaspect([1 1 1])
plot([.1 .1],ylim,'r--')
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 4 4];
fig.PaperPositionMode = 'Manual';
fig.PaperSize = [4 4];
if savefig
    fn = ['distribution_switch_times_' which_switch correction_str '.svg'];
    print(fig, fullfile(dp.sta_fig_dir, fn),'-dsvg')
    fn = ['switch_times_' which_switch correction_str '.mat'];
    save(fullfile(dp.sta_fig_dir, fn), 'middle_ts');
end

figure(6); clf; hold on;
%plot(sort(middle_ts),1:length(middle_ts),'k-','linewidth',2)
histogram(middle_ts,35,'facecolor','k')
plot([0 0], ylim, 'k-','linewidth',2)
ylabel('Count')
xlabel(['Switch time (s)'])
set(gca,'fontsize',16)
pbaspect([1 1 1])
plot([.1 .1],ylim,'r-','linewidth',2)
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 4 4];
fig.PaperPositionMode = 'Manual';
fig.PaperSize = [4 4];
if savefig
    fn = ['distribution_switch_times_histogram_' which_switch correction_str '.svg'];
    print(fig, fullfile(dp.sta_fig_dir, fn),'-dsvg')
end


%%
plot_val = (1-ptiles).*(time_vec<=0) + ptiles.*(time_vec>0);
%plot_val = 1-ptiles;
a = .05;
nconsec = 2


sigp = ptiles < a | ptiles > (1-a);
csigp = cumsum(sigp,2);

postswitchsigp = sigp.*(time_vec>0);
movsumsig = movsum(postswitchsigp,nconsec,2);
consecsig = movsumsig==nconsec;
close(figure(10));
fh = figure(10); clf
set(fh,'position',[2 2 6 4],'papersize',[ 6 4])

switch 1
    case 0
        consecsig = consecsig(good_ind,:);
        %consecsig(:,end) = nconsec;
        plot_val = plot_val(ind,:);
    case 1
        consecsig = consecsig(good_ind,:);
        consecsig(:,end) = nconsec;
        plot_val = plot_val(good_ind,:);
    case 2
end
[~, i_s, j_s] = sort_by_peak(consecsig(:,:));

plotMat = plot_val(i_s(j_s>1),:);


ax = axes


n = 20;
sig_bins = [linspace(0, .05,n) linspace(.1,.9,n) linspace(.95, 1,n)];
sig_val = sig_bins(1:end-1) + diff(sig_bins)/2;
cm = flipud(colormapLinear(dp.left_color,n-1));
cm = [cm; ones(n,3);  colormapLinear(dp.right_color,n-1)];

[M, i, plotBin] = histcounts(plotMat,sig_bins);
%plotBin = sig_val(plotBin);
plotMat = reshape(plotBin,size(plotMat,1),size(plotMat,2));
imagesc(plotMat,'CDataMapping','direct','x',time_vec)

colormap(cm)
%colormap(jet)

caxis('auto')
hold on
plot([0 0],ylim,'k--')
cb = colorbar('eastoutside')
xlabel('time from state switch (s)')
drawnow
axpos = get(ax,'position');

set(cb, 'position', cb.Position + [.07 0 0 -.4]);
set(ax,'position',axpos+ [.00 0 .065 0])
drawnow
cb.YTick = [n-.5 size(cm,1)-n+.5];
cb.YTickLabel = [.05 .95];
ax.TickDir = 'out';
title(cb,{'d'' %tile'})
box(ax,'off')

ex_y = sort([find(cellids(i_s(j_s>1)) == 18181) ...
    find(cellids(i_s(j_s>1)) == 17784) ...
    find(cellids(i_s(j_s>1)) == 16857)]);
%ax.YTick = ex_y;
%ax.YTickLabel = [];

ylabel('cells with 10 significant conecutive timepoints')
ylabel({'cell #' '(sorted by selectivity' 'onset post-switch)'})

xlim([-.45 .75])

print(fh, fullfile(dp.fig_dir, 'popSTA'), '-dsvg', '-painters')

%%
figure(11); clf

imagesc(consecsig(i_s(j_s>1),:),'x',time_vec)
i_s = i_s ;

hold on
plot([0 0],ylim,'k--')
colorbar
title('p<.05 for 5 consecutive timepoints','fontsize',10,'fontweight','normal')
% 
% 
% plot(0, ex_y,'o','markersize',5,'markerfacecolor','m','markeredgecolor','k')
% 
% ex_y = find(cellids(i_s(j_s>1)) == 16857);
% plot(0, ex_y,'o','markersize',5,'markerfacecolor','c','markeredgecolor','k')
% 
% ex_y = find(cellids(i_s(j_s>1)) == 17784);
% plot(0, ex_y,'o','markersize',5,'markerfacecolor','g','markeredgecolor','k')

xlim([-.5 .75])


keyboard
%%
subplot(212)
diffSTA   = leftSTA - rightSTA;
imagesc(diffSTA(i_s(j_s>1),:),'x',res{1}.lags)

%%

diffSTA   = leftSTA - rightSTA;
sig_level = sum(sigp,2);

which_cells = 1:ncells;
% [~, which_cells] = sort(sig_level);
% which_cells = find(sig_level > 20)

%which_cells = [21  510 474];
fh = figure(10); clf
set(fh, 'position',[5 5 3 6],'papersize',[3 3])
ax = subplot(211);
diffSTA = diffSTA(i_s(j_s>1),:);
NA = ~isnan(diffSTA);
%diffSTA = sort_by_peak(diffSTA,-abs(diffSTA))
imagesc(diffSTA,'x',res{1}.lags,'Alphadata',NA,'parent',ax)
set(ax,'color',[1 1 1].*.8)

caxis([-1 1].*5)
cm = flipud(colormapLinear(dp.left_color,50));
cm = [cm; colormapLinear(dp.right_color,50)];
colormap(colormapRedBlue)
colormap(cm)
xlim([-.5 .75])
hold on
plot([0 0],ylim,'k--')
%%

subplot(212)
imagesc(ptiles(which_cells,:),'x',time_vec)
caxis([0 1])
%%
% some debugging plotting code
if 0
dex = 1;
computed_cells = cellids(logical(computed));
%% Double check consistency of indexes and such
plotSTA(computed_cells(dex));
drawnow;
pause(1)
plotDprime(computed_cells(dex));
hold on;
% dprime values match
plot(time_vec, dprime(dex,:),'m')
% pvalues match
pvals_match = sum(res{dex}.pval(pval_plot_lags)' == pvals(dex,:));
plot(time_vec, pvals(dex,:),'m')
ylim([-1 1.5])

% pvals indices match
 pplot = pvals(dex,:);
 pplot(pplot < 0.95 & pplot > 0.05) = NaN;
 pplot1 = pplot;
 pplot1(pplot < 0.95) = NaN;
 pplot1(~isnan(pplot1)) = 1.2;
 pplot2 = pplot;
 pplot2(pplot > 0.05) = NaN;
 pplot2(~isnan(pplot2)) = 1.2;
 plot(res{dex}.lags(pval_plot_lags),pplot1, 'gs', 'linewidth',1, 'markerfacecolor','g','markersize',5);
 plot(res{dex}.lags(pval_plot_lags),pplot1, 'g-', 'linewidth',5);
 plot(res{dex}.lags(pval_plot_lags),pplot2, 'bs', 'linewidth',1, 'markerfacecolor','b','markersize',5);
 plot(res{dex}.lags(pval_plot_lags),pplot2, 'b-', 'linewidth',5);

% plot map indicies match
 pmplot = plot_map(dex,:);
 pmplot(pmplot < 0.95 & pmplot > 0.05) = NaN;
 pmplot1 = pmplot;
 pmplot1(pmplot < 0.95) = NaN;
 pmplot1(~isnan(pmplot1)) = 1;
 pmplot2 = pmplot;
 pmplot2(pmplot > 0.05) = NaN;
 pmplot2(~isnan(pmplot2)) = 1;
 plot(res{dex}.lags(pval_plot_lags),pmplot1, 'gs', 'linewidth',1, 'markerfacecolor','g','markersize',5);
 plot(res{dex}.lags(pval_plot_lags),pmplot1, 'g-', 'linewidth',5);
 plot(res{dex}.lags(pval_plot_lags),pmplot2, 'bs', 'linewidth',1, 'markerfacecolor','b','markersize',5);
 plot(res{dex}.lags(pval_plot_lags),pmplot2, 'b-', 'linewidth',5);

plot(time_vec, Q1(dex,:).*.9,'k--')
plot(time_vec, Q2(dex,:).*.8,'k--')
plot(time_vec, Q3(dex,:).*.7,'k--')
plot(time_vec, Q4(dex,:).*.6,'k--')

end

