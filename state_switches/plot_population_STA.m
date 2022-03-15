function [axsta, pop_res ] = plot_population_STA(varargin)
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
p               = p.Results;
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
cell_list   = dyn_cells_db;
select_str  = 'normmean > 0' ;

psr         = cell2mat(extracting(cell_list, 'prefsideright', select_str));
normmean    = cell2mat(extracting(cell_list, 'normmean', select_str));

leftstd     = nanstd(leftSTA,[],2);
rightstd    = nanstd(rightSTA,[],2);

leftdenom   = leftstd;
rightdenom  = rightstd;

leftSTAz    = leftSTA ./ leftdenom;
rightSTAz   = rightSTA ./ rightdenom;

good_cells  = cf.good_cells;

psr = psr == 1;
STR_pref = leftSTAz;
STR_pref(psr,:) = rightSTAz(psr ,:);
STR_npref  = rightSTAz;
STR_npref(psr,:) = leftSTAz(psr ,:);

STR_pref_good = STR_pref(good_cells,:);
STR_npref_good = STR_npref(good_cells,:);
%% plot population average traces
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

STR_npref_mn = nanmean(STR_npref_good);
STR_npref_sem = nansem(STR_npref_good);
STR_pref_mn = nanmean(STR_pref_good);
STR_pref_sem = nansem(STR_pref_good);

axsta.TickDir = 'out';
% plot switch to non-preferred side
shadedErrorBar(plot_lags(negdex),STR_npref_mn(:,negdex),...
    STR_npref_sem(:,negdex),{'color', dp.pref_color,'parent',axsta},[],sat)
shadedErrorBar(plot_lags(posdex),STR_npref_mn(:,posdex),...
    STR_npref_sem(:,posdex),{'color', dp.npref_color},[],sat-satdiff)
% plot switch to preferred side
shadedErrorBar(plot_lags(negdex),STR_pref_mn(:,negdex),...
    STR_pref_sem(:,negdex),{'color', dp.npref_color},[],sat)
shadedErrorBar(plot_lags(posdex),STR_pref_mn(:,posdex),...
    STR_pref_sem(:,posdex),{'color', dp.pref_color},[],sat-satdiff)

axis tight
ylim([-1 1]*.425)
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
set(fh,'position',[3 11 ppos],'papersize',[ppos],'paperpositionmode','auto');

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
axpos = get(gca,'position');
cb = colorbar;
set(gca,'position',axpos);
cb.Position = cb.Position + [-.03 0 -.01 -.4];
title('selective cells','fontweight','normal')
figsavefname = ['population_' which_switch '_STA_heat_pref'];
print(fh, fullfile(dp.fig_dir, figsavefname) ,fig_type,'-painters')
%% heat map side - not used
fh= figure(3); clf
fht = 3;
fw  = 3.75;
ppos = [fw fht];
set(fh,'position',[3 7 ppos],'papersize',[ppos],'paperpositionmode','auto')

%STR_diff = STR_pref_good - STR_npref_good;
STR_diffz = rightSTAz - leftSTAz;
STR_diffz = STR_diffz(good_cells,:);
STR_diffz = [STR_diffz(psr(good_cells),:); STR_diffz(~psr(good_cells),:)];
STR_diffz(:,plot_lags < 0) = -STR_diffz(:,plot_lags < 0) ;
imagesc(STR_diffz,'x',plot_lags)


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
axpos = get(gca,'position');
nc = 10;
mid_color = [1 1 1].*1;
cm = flipud(colormapLinear(dp.right_color,nc,mid_color));
cm2 = colormapLinear(dp.left_color,nc,mid_color);

clrs = [cm(1:end-1,:); cm2(1:end,:)];
colormap(flipud(clrs))
%colormap(colormapRedBlue)
cb = colorbar;
set(gca,'position',axpos);
cb.Position = cb.Position + [-.03 0 -.01 -.4];
figsavefname = ['population_' which_switch '_STA_heat_side'];
print(fh, fullfile(dp.fig_dir, figsavefname) ,fig_type,'-painters')
%% heat map split into two plots by preference - not used
fh= figure(3); clf
fht = 2;
fw  = 3.75;
ppos = [fw fht];
set(fh,'position',[3 7 ppos],'papersize',[ppos],'paperpositionmode','auto')

STR_diffz = rightSTAz - leftSTAz;
STR_diffz = STR_diffz(good_cells,:);
STR_diffz(:,plot_lags < 0) = -STR_diffz(:,plot_lags < 0) ;

ax1 = subplot(211);
imagesc(STR_diffz(psr(good_cells),:),'x',plot_lags)
caxis([-1 1]*1)
xlim([ min_t max_t])
ax2 = subplot(212);
imagesc(STR_diffz(~psr(good_cells),:),'x',plot_lags)
caxis([-1 1]*1)
xlim([ min_t max_t])
nr = sum(psr(good_cells));
nl = sum(~psr(good_cells));

ylabel(ax1,'right cell #')
ylabel(ax2,'left cell #')
xlabel(ax2,['time from ' which_switch ' state change (s)'])
drawnow

ax1pos = get(ax1,'position');
ax2pos = get(ax2,'position');
set(ax2,'position',ax2pos+[0 0 0 .1]);
ax2pos = get(ax2,'position');

ax1pos = get(ax1,'position');
set(ax1,'position',[ax1pos([1 2]) ax2pos(3) nr/nl*ax2pos(4)]-[0 .05 0 0]);
%set(ax1,'position',[ax2pos([1 2]) ax1pos(3) nl/nr*ax1pos(4)])
set(ax1,'xticklabel',[],'ytick',[1 30 60]);
set(ax2,'xticklabel',[],'ytick',[1 30 60]);
box(ax1,'off')
box(ax2,'off')

title(ax1,'STR for all selective cells','fontweight','normal')
figsavefname = ['population_' which_switch '_STA_heat_side_split'];
print(fh, fullfile(dp.fig_dir, figsavefname) ,fig_type,'-painters')

%%
% make map of time points for which each cell is encoding signficantly

pval_plot_lags  = plot_lags > min_t & plot_lags < max_t;
ptiles_plot     = ptiles(:,pval_plot_lags);
posdex          = posdex(pval_plot_lags);
negdex          = negdex(pval_plot_lags);
alpha           = p.alpha;

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
plot_map(isnan(plot_map)) = 0.5;
temp = plot_map - 0.5;
temp(temp < threshold & temp > -threshold) = 0;
plot_map = temp + 0.5;
sort_dex = sum(temp,2);
[~, sort_ind] = sort(sort_dex,1,'descend');
temp                = abs(temp(good_ind,:));
bin_temp            = temp;
bin_temp(temp > 0)  = 1;

% look at fraction of significant cells

average_siggy   = sum(bin_temp,1)./size(temp,1);
n               = size(temp,1);
p               = movmean(average_siggy,3);
se              = 1.96.*sqrt((p.*(1-p))./n);
time_vec        = plot_lags(pval_plot_lags);

datasavefname   = ['fraction_data_' which_switch correction_str '.mat'];
% Save data from this analysis
save(fullfile(dp.sta_dir, datasavefname), 'time_vec','average_siggy','n')
%%
pop_res.STR_time = plot_lags;
pop_res.STR_npref_mn = STR_npref_mn;
pop_res.STR_npref_sem = STR_npref_sem;
pop_res.STR_pref_mn  = STR_pref_mn;
pop_res.STR_pref_sem = STR_pref_sem;
pop_res.STR_diff = STR_diff;
