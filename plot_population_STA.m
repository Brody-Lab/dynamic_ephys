function plot_population_STA(which_switch, savefig,correction, bad_strength) 


%which_switch = 'generative';
if nargin < 2
    savefig = 0;
end
if nargin < 3
    correction = 2;
end
if nargin < 4
    bad_strength = 0;
end
[res, dprime, pvals,cellids, computed,not_computed] = load_STA_population(which_switch, 0,1, bad_strength);

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
    m = numel(pvals);
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
plot_map = pvals;
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
    print(fig, ['../figures/STA/' which_switch 'switch_encoding' correction_str],'-depsc','-opengl')
end

%%%%%%%
if strcmp(which_switch,'model')
    load('ephys_summary/sorted_cells.mat')
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
        print(fig, ['../figures/STA/' which_switch 'switch_encoding' correction_str '_sequence_labels'],'-depsc','-opengl')
    end
end


% look at fraction of significant cells
temp                = abs(temp);
bin_temp            = temp;
bin_temp(temp > 0)  = 1;
average_siggy       = sum(bin_temp,1)./size(temp,1);
n = size(temp,1);
p = movmean(average_siggy,3);
se = 1.96.*sqrt((p.*(1-p))./n);
time_vec            = res{1}.lags(pval_plot_lags);
figure(2); clf;
shadedErrorBar(time_vec, average_siggy.*100, se.*100, 'k')
%plot(time_vec, movmean(average_siggy.*100,3), 'k','linewidth',2)
ylim([0 20])
hold on
plot([.1 .1], ylim, 'r--')
plot([0 0], [-1 100], 'k--')
set(gca, 'fontsize',16)
ylabel('Significant cells (% of total)')
xlabel(['Time from ' which_switch ' switch (s)'])
pbaspect([1 1 1]);
xlim([-.41 .41])
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 4 4];
fig.PaperPositionMode = 'Manual';
fig.PaperSize = [4 4];

if savefig
print(fig, ['../figures/STA/fraction_' which_switch 'switch_encoding' correction_str '.svg'],'-dsvg')

% Save data from this analysis
save(['../figures/STA/fraction_data_' which_switch correction_str '.mat'], 'time_vec','average_siggy','n')
end

%%%% look at distribution of dprime values
figure(3); clf;
dabs = abs(dprime);
mean_dabs = nanmean(dabs,2);
histogram(mean_dabs(mean_dabs < 10),50, 'facecolor','k')
ylim([0 inf])%
xlim([0 .5])
ylabel('count')
xlabel('mean d'' per cell')
set(gca,'fontsize',16)
pbaspect([1 1 1])
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 4 4];
fig.PaperPositionMode = 'Manual';
fig.PaperSize = [4 4];
if savefig
    print(fig, ['../figures/STA/distribution_dprime_' which_switch '.svg'],'-dsvg')
end

%%%% look for consistent encoding cells
const_map   = pvals;
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
    print(fig, ['../figures/STA/' which_switch '_switch_labeled_switch_encoding' correction_str],'-depsc','-opengl')
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
    print(fig, ['../figures/STA/distribution_switch_times_' which_switch correction_str '.svg'],'-dsvg')
    save(['../figures/STA/switch_times_' which_switch correction_str '.mat'], 'middle_ts');
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
    print(fig, ['../figures/STA/distribution_switch_times_histogram_' which_switch correction_str '.svg'],'-dsvg')
end

keyboard
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

