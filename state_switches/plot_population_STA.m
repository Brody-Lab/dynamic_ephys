function [fh ] = plot_population_STA(varargin)
p = inputParser;
addParameter(p, 'which_switch', 'model')
addParameter(p, 'savefig', 1)
addParameter(p, 'correction_num', 2)
addParameter(p, 'recompute', 0)
addParameter(p, 'lag', 0)
addParameter(p, 'max_t', 1)
addParameter(p, 'min_t', -1)
addParameter(p, 'fig_type', '-dsvg')
addParameter(p, 'clear_bad_strengths',1 )
addParameter(p, 'bad_strength', 0)
addParameter(p, 't_buffers', [0 0])
addParameter(p, 'min_post_dur', 0)
addParameter(p, 'min_pre_dur', 0)
addParameter(p, 'alpha', 0.05)
parse(p,varargin{:})
p = p.Results;
which_switch    = p.which_switch;
savefig         = p.savefig;
correction      = p.correction_num;
bad_strength    = p.bad_strength;
force           = p.recompute;
lag             = p.lag;
dp              = set_dyn_path;
max_t           = p.max_t;
min_t           = p.min_t;
fig_type        = p.fig_type;
[res, cellids, computed, not_computed] = ...
    load_STA_population(which_switch, 'force',force,'slim_data',1,...
    'bad_strength', p.bad_strength, ...
    't_buffers', p.t_buffers, ...
    'min_pre_dur', p.min_pre_dur,'min_post_dur', p.min_post_dur,...
    'clear_bad_strengths',p.clear_bad_strengths);
dprime = [];
ptiles = [];

for ccall = 1:length(cellids)
    if ~isempty(res{ccall})
        pval_plot_lags = res{ccall}.lags > -Inf & res{ccall}.lags < Inf;
        dprime = [dprime; res{ccall}.dprime_real(pval_plot_lags)'];
        ptiles =  [ptiles;  res{ccall}.pval(pval_plot_lags)'];
        cc = ccall;
    end
end
%%
ncells      = length(cellids);
nlags       = length(res{1}.lags);
plot_lags   = res{1}.lags - p.lag;
leftSTA     = nan(ncells,nlags);
rightSTA    = nan(ncells,nlags);
for cc = 1:length(cellids)
    leftSTA(cc,:)   = res{cc}.STR_left_real;% / nanstd(res{cc}.STR_left_real);
    rightSTA(cc,:)  = res{cc}.STR_right_real;% / nanstd(res{cc}.STR_right_real); 
end
%%
cout_auc_file = fullfile(dp.ephys_summary_dir,'cout_auc.mat');
cf = load(cout_auc_file,'cellids','good_cells');
good_ind = ismember(cellids,cf.cellids(cf.good_cells));

%% population STA traces
cell_list = dyn_cells_db;
select_str = 'normmean > 0' ;

psr = cell2mat(extracting(cell_list, 'prefsideright', select_str))
normmean = cell2mat(extracting(cell_list, 'normmean', select_str));

leftstd = nanstd(leftSTA,[],2);
rightstd = nanstd(rightSTA,[],2);

switch 1
    case 0
        leftdenom = normmean;
        rightdenom = normmean;
    case 1
        leftdenom = leftstd;
        rightdenom = rightstd;
end

leftSTAz = leftSTA ./ leftdenom;
rightSTAz = rightSTA ./ rightdenom;

good_cells = true(size(cf.good_cells));
good_cells = cf.good_cells;

psr = psr == 1;
STR_pref = leftSTAz;
STR_pref(psr,:) = rightSTAz(psr ,:);
STR_npref  = rightSTAz;
STR_npref(psr,:) = leftSTAz(psr ,:);

STR_pref_good = STR_pref(good_cells,:);
STR_npref_good = STR_npref(good_cells,:);
%%
fh = figure(1); clf
fht = 2;
fw  = 3.75;
ppos = [fw fht];
set(fh,'position',[2 5 fw fht],'papersize', [fw fht])

axsta = axes;
hold(axsta,'on')

posdex = plot_lags > -0.001 & plot_lags < max_t;
negdex = plot_lags <  0.001 & plot_lags > min_t;
sat = .7; % .9
satdiff = 0; % .2

axsta.TickDir = 'out';
shadedErrorBar(plot_lags(negdex),nanmean(STR_npref_good(:,negdex)),...
    nansem(STR_npref_good(:,negdex)),{'color', dp.pref_color,'parent',axsta},[],sat)
shadedErrorBar(plot_lags(negdex),nanmean(STR_pref_good(:,negdex)),...
    nansem(STR_pref_good(:,negdex)),{'color', dp.npref_color},[],sat)
shadedErrorBar(plot_lags(posdex),nanmean(STR_npref_good(:,posdex)),...
    nansem(STR_npref_good(:,posdex)),{'color', dp.npref_color},[],sat-satdiff)
shadedErrorBar(plot_lags(posdex),nanmean(STR_pref_good(:,posdex)),...
    nansem(STR_pref_good(:,posdex)),{'color', dp.pref_color},[],sat-satdiff)
axis tight
ylim([-1 1]*.4)
set(axsta, 'ytick', [-.4 -.2 0 .2 .4])
plot([0 0],ylim,'--k')
xlim([min_t max_t])


%cb = colorbar
xlabel(['time from ' which_switch ' state change (s)'])
ylabel('\Delta FR (z-score)')
title('selective population average STR','fontweight','normal')
drawnow

figsavefname = ['population_' which_switch '_STA'];
print(fh, fullfile(dp.fig_dir, figsavefname) ,fig_type,'-painters')
%% heat map pref
fh= figure(2); clf
set(fh,'position',[3 11 ppos],'papersize',[ppos],'paperpositionmode','auto')

STR_diff = STR_pref_good - STR_npref_good;
%STR_diff = STR_pref - STR_npref;
imagesc((STR_diff),'x',plot_lags)
colormap(colormapRedBlue)
caxis([-1 1]*1.5)

xlim([ min_t max_t])
hold on
plot([0 0],ylim,'k')
ylabel('cell #')
xlabel(['time from ' which_switch ' state change (s)'])
box off
drawnow
axpos = get(gca,'position')
cb = colorbar
set(gca,'position',axpos)
cb.Position = cb.Position + [-.03 0 -.01 -.4]
title('selective cells','fontweight','normal')
figsavefname = ['population_' which_switch '_STA_heat_pref'];
print(fh, fullfile(dp.fig_dir, figsavefname) ,fig_type,'-painters')
%% heat map side
fh= figure(3); clf
fht = 3;
fw  = 3.75;
ppos = [fw fht];
set(fh,'position',[3 7 ppos],'papersize',[ppos],'paperpositionmode','auto')

%STR_diff = STR_pref_good - STR_npref_good;
STR_diff = rightSTAz - leftSTAz;
STR_diff = STR_diff(good_cells,:)
STR_diff = [STR_diff(psr(good_cells),:); STR_diff(~psr(good_cells),:)];
STR_diff(:,plot_lags < 0) = -STR_diff(:,plot_lags < 0) ;
imagesc(STR_diff,'x',plot_lags)


hold on
plot( xlim, [1 1]*sum(psr(good_cells)),'k')
caxis([-1 1]*1)
xlim([ min_t max_t])
hold on
plot([0 0],ylim,'k')
ylabel('cell #')
xlabel(['time from ' which_switch ' state change (s)'])
box off
title('STR for all selective cells','fontweight','normal')
drawnow
axpos = get(gca,'position')
nc = 10
mid_color = [1 1 1].*1;
cm = flipud(colormapLinear(dp.right_color,nc,mid_color));
cm2 = colormapLinear(dp.left_color,nc,mid_color);

clrs = [cm(1:end-1,:); cm2(1:end,:)];
colormap(flipud(clrs))
%colormap(colormapRedBlue)
cb = colorbar
set(gca,'position',axpos)
cb.Position = cb.Position + [-.03 0 -.01 -.4]
figsavefname = ['population_' which_switch '_STA_heat_side'];
print(fh, fullfile(dp.fig_dir, figsavefname) ,fig_type,'-painters')
%%
fh= figure(3); clf
fht = 2;
fw  = 3.75;
ppos = [fw fht];
set(fh,'position',[3 7 ppos],'papersize',[ppos],'paperpositionmode','auto')

%STR_diff = STR_pref_good - STR_npref_good;
STR_diff = rightSTAz - leftSTAz;
STR_diff = STR_diff(good_cells,:)
STR_diff(:,plot_lags < 0) = -STR_diff(:,plot_lags < 0) ;


ax1 = subplot(211)
imagesc(STR_diff(psr(good_cells),:),'x',plot_lags)
caxis([-1 1]*1)
xlim([ min_t max_t])
ax2 = subplot(212)
imagesc(STR_diff(~psr(good_cells),:),'x',plot_lags)
caxis([-1 1]*1)
xlim([ min_t max_t])
nr = sum(psr(good_cells));
nl = sum(~psr(good_cells));


ylabel(ax1,'right cell #')
ylabel(ax2,'left cell #')
xlabel(ax2,['time from ' which_switch ' state change (s)'])
drawnow

ax1pos = get(ax1,'position')
ax2pos = get(ax2,'position')
set(ax2,'position',ax2pos+[0 0 0 .1])
ax2pos = get(ax2,'position')

ax1pos = get(ax1,'position')
set(ax1,'position',[ax1pos([1 2]) ax2pos(3) nr/nl*ax2pos(4)]-[0 .05 0 0])
%set(ax1,'position',[ax2pos([1 2]) ax1pos(3) nl/nr*ax1pos(4)])
set(ax1,'xticklabel',[],'ytick',[1 30 60])
set(ax2,'xticklabel',[],'ytick',[1 30 60])
box(ax1,'off')
box(ax2,'off')

title(ax1,'STR for all selective cells','fontweight','normal')
figsavefname = ['population_' which_switch '_STA_heat_side_split'];
print(fh, fullfile(dp.fig_dir, figsavefname) ,fig_type,'-painters')

%%

% make map of time points for which each cell is encoding signficantly
%figure(1); clf
figure(121); clf
pval_plot_lags = plot_lags > min_t & plot_lags < max_t;
ptiles_plot = ptiles(:,pval_plot_lags)
clrs = color_set(6);
posdex = (plot_lags > -0.0001);
negdex = (plot_lags <  0.0001);
posdex = posdex(pval_plot_lags);
negdex = negdex(pval_plot_lags);
cmap = [repmat(clrs(1,:),10,1); repmat([1 1 1], 200, 1); repmat(clrs(end,:),10,1)];

alpha = p.alpha;

if correction == 1 % bonferroni corrections
    m = numel(ptiles_plot);
    threshold = 0.5 - alpha/2/ m;
    correction_str = '_bonferroni';
elseif correction == 2 % modified bonferroni, because tests are correlated
    % threshold determined by emperically checking for false positives...
    %in the shuffle data, and asking what index keeps it below 0.05
    m = 5;
    threshold = 0.5 - alpha/2/ m;
    correction_str = '_bonferroni_modified';
else
    threshold = 0.5 - alpha/2;
    correction_str = '';
end

correction_str = [correction_str '_' num2str(bad_strength)];
plot_map = ptiles_plot;
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
    imagesc(sorted_plot_map,'x',plot_lags(pval_plot_lags),[0 1]); 
else
    imagesc(plot_map(sort_ind,:),'x',plot_lags(pval_plot_lags),[0 1]); 
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
%average_siggy = average_siggy;
n = size(temp,1);
p = movmean(average_siggy,3);
se = 1.96.*sqrt((p.*(1-p))./n);
time_vec            = plot_lags(pval_plot_lags);
%%
figure(3); clf;
ax = axes;
shadedErrorBar(time_vec, average_siggy.*100, se.*100, 'k')
%plot(time_vec, movmean(average_siggy.*100,3), 'k','linewidth',2)
ylim([0 50])
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
    figsavefname = ['fraction_' which_switch 'switch_encoding' correction_str ];
    print(fig, fullfile(dp.fig_dir, figsavefname) ,fig_type,'-painters')
    datasavefname = ['fraction_data_' which_switch correction_str '.mat'];
    % Save data from this analysis
    save(fullfile(dp.sta_fig_dir, datasavefname), 'time_vec','average_siggy','n')
end
%%
if 0
%%%% look at distribution of dprime values
figure(4); clf;
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
    figsavefname = ['distribution_dprime_' which_switch ];
    print(fig, fullfile(dp.sta_fig_dir, figsavefname),fig_type)
end
%%

%%%% look for consistent encoding cells
const_map   = ptiles_plot(:,:);
const_map(const_map > 0.99)     = 1;
const_map(const_map < 0.01)     = -1;
const_map(const_map < 1 & const_map > -1) = 0;
const_pmap   = plot_map;
const_pmap(const_pmap > (1-alpha/2))     = 1;
const_pmap(const_pmap < (alpha/2))     = -1;
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
plot_map(plot_map > (1-alpha/2)) =1;
plot_map(plot_map < alpha/2) =0;
if sortplot
    imagesc(plot_map(sort_ind,:),'x',plot_lags(pval_plot_lags),[0 1]); 
%    sorted_plot_map = plot_map(sort_ind,:);
%    all_zero_dex = all(sorted_plot_map == 0.5,2);
%    sorted_plot_map(all_zero_dex,:) = [];
%    imagesc(sorted_plot_map,'x',plot_lags(pval_plot_lags),[0 1]); 
else
    imagesc(plot_map,'x', plot_lags(pval_plot_lags),[0 1]); colormap(cmap);
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
    fn = ['distribution_switch_times_' which_switch correction_str ];
    print(fig, fullfile(dp.sta_fig_dir, fn),fig_type)
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
    fn = ['distribution_switch_times_histogram_' which_switch correction_str ];
    print(fig, fullfile(dp.sta_fig_dir, fn),fig_type)
end




%%
selective_only = 1;
plot_val = (1-ptiles_plot).*(time_vec<=0) + ptiles_plot.*(time_vec>0);
%plot_val = 1-ptiles;
a = .05;

nconsec = 40;
use_pre  = 0;
use_mid = 1

sigp = ptiles_plot < alpha/2 | ptiles_plot > (1-alpha);
csigp = cumsum(sigp,2);

postswitchsigp = sigp.*(time_vec>0);
preswitchsigp = sigp.*(time_vec<0.50);
    
postswitchsigp = sigp.*(time_vec>0);
postmovsumsig = movsum(postswitchsigp,nconsec,2);
postconsecsig = postmovsumsig==nconsec;

premovsumsig = movsum(preswitchsigp,nconsec,2);
preconsecsig = cumsum(premovsumsig>=nconsec,2);
   
ylab = 'cell #';
    
if use_pre
    consecsig = preconsecsig;
    movsumsig = premovsumsig;
    ylab = {'cell #' 'sorted pre-switch selectivity offset'}
elseif use_mid
    consecsig = preconsecsig + postconsecsig;
    movsumsig = premovsumsig + postmovsumsig;
    consecsig = (movsum(~sigp(:,time_vec > -.2 & time_vec < .1),3,2));
    consecsig = (movsum(ptiles_plot(:,time_vec > -.2 & time_vec < .1)-.5,15,2));
    %consecsig = (nansum(ptiles_plot(:,time_vec > -.2 & time_vec < .1),2));
else
    consecsig = postconsecsig;
    movsumsig = postmovsumsig;
    ylab = {'cell #' 'sorted post-switch selectivity onset'};

end
close(figure(10));
fh = figure(10); clf
set(fh,'position',[2 2 6 5],'papersize',[ 6 4])
ax = axes
set(ax,'position',[.15 .12 .75 .8])
if ~selective_only
    plot_ind = true(size(consecsig,1),1);
else
    plot_ind = good_ind;
end

%consecsig = consecsig(plot_ind,:);

switch 1
    case 0
        consecsig = consecsig(plot_ind,:);
        %consecsig(:,end) = nconsec;
        plot_val = plot_val(plot_ind,:);
        min_val = 1;
    case 1
        consecsig = consecsig(plot_ind,:);
        min_val = -Inf;
        plot_val = plot_val(plot_ind,:);
%     case 3
%         consecsig = cumsum(movsumsig>=nconsec,2);
% 
%         consecsig = cumsum(postswitchsigp,2)
%     case 4
%         consecsig = cumsum(movsumsig>=nconsec,2);

                
end
[~, i_s, j_s] = sort_by_peak(consecsig(:,:));

plotMat = plot_val(i_s(j_s>min_val),:);

n = 20;
sig_bins = unique([linspace(0, alpha/2,n) ...
    linspace(alpha/2,1-alpha/2,n) ...
    linspace(1-alpha/2, 1,n)]);
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
xlabel(['time from ' which_switch ' switch (s)'])
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

ylabel(ylab)


print(fh, fullfile(dp.fig_dir, [which_switch 'popSTA']), fig_type, '-painters')

%%
diffSTA   = leftSTA - rightSTA;
diffSTA = -diffSTA .* (2*(plot_lags > 0)-1)

fh = figure(100); clf
set(fh, 'position',[5 5 4 4],'papersize',[3 3])
which_cells = i_s(j_s>1);
which_cells = 1:size(diffSTA,1)
if selective_only
    which_cells = good_ind;
end
pre_mn = nanmean(diffSTA(which_cells,plot_lags<0 & plot_lags>-.4),2);
post_mn = nanmean(diffSTA(which_cells,plot_lags>0& plot_lags<.5),2);
r_cells = pre_mn > 0 & post_mn > 0;
l_cells = pre_mn < 0 & post_mn < 0;
b_cells = ~r_cells & ~l_cells;
plot(pre_mn(r_cells),post_mn(r_cells),'.',...
    'markersize',10,'color',dp.right_color)
hold on
plot(pre_mn(l_cells),post_mn(l_cells),'.',...
    'markersize',10,'color',dp.left_color)
plot(pre_mn(b_cells),post_mn(b_cells),'.',...
    'markersize',10,'color',[.9 .6 .9])
axis image
r = max(abs([pre_mn; post_mn]));
xlim([-1 1]*r)
ylim([-1 1]*r)
plot([0 0],ylim,'k')
plot(xlim,[0 0],'k')

box off
xlabel('pre-switch rate difference')
ylabel('post-switch rate difference')
axis('tight')
fn = ['rate_diff_scatter' which_switch correction_str ];

print(fh, fullfile(dp.sta_fig_dir, fn),fig_type,'-painters')
%%
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
 pplot(pplot < (1-alpha/2) & pplot > alpha/2) = NaN;
 pplot1 = pplot;
 pplot1(pplot < (1-alpha/2)) = NaN;
 pplot1(~isnan(pplot1)) = 1.2;
 pplot2 = pplot;
 pplot2(pplot > alpha/2) = NaN;
 pplot2(~isnan(pplot2)) = 1.2;
 plot(res{dex}.lags(pval_plot_lags),pplot1, 'gs', 'linewidth',1, 'markerfacecolor','g','markersize',5);
 plot(res{dex}.lags(pval_plot_lags),pplot1, 'g-', 'linewidth',5);
 plot(res{dex}.lags(pval_plot_lags),pplot2, 'bs', 'linewidth',1, 'markerfacecolor','b','markersize',5);
 plot(res{dex}.lags(pval_plot_lags),pplot2, 'b-', 'linewidth',5);

% plot map indicies match
 pmplot = plot_map(dex,:);
 pmplot(pmplot < (1-alpha/2) & pmplot > alpha/2) = NaN;
 pmplot1 = pmplot;
 pmplot1(pmplot < (1-alpha/2)) = NaN;
 pmplot1(~isnan(pmplot1)) = 1;
 pmplot2 = pmplot;
 pmplot2(pmplot > alpha/2) = NaN;
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
end
