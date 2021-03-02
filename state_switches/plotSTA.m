function res = plotSTA(cellid, varargin)
p = inputParser;
addParameter(p,'which_switch', 'model')
addParameter(p,'save_fig', 1)
addParameter(p,'fig_num', [])
addParameter(p,'plot_residual',1)
addParameter(p,'lag',0)
addParameter(p,'max_t',.75)
addParameter(p,'min_t',-.75)
addParameter(p,'fig_type','svg')


parse(p,varargin{:})
p = p.Results;
which_switch    = p.which_switch;
save_fig        = p.save_fig;
max_t = p.max_t;
min_t = p.min_t;
fig_type = ['-d' p.fig_type];

if isempty(p.fig_num)
    fh = figure; clf;
    fpos = [4 4]
else
    fh = figure(p.fig_num); clf
    fpos = get(fh,'position')
end
set(fh,'position',[fpos([1 2]) 6 3],'papersize',[6 3],'paperpositionmode','auto')

dp = set_dyn_path;
if isstruct(cellid)
    res = cellid;
    cellid = res.cellid;
    which_switch = res.which_switch;
else
    res =  compute_switch_triggered_average(cellid,'post',2,...
        'which_switch',which_switch, 'n_shuffles', 1000,...
        'save_file',1,'mask_other_switch',1,'force',1);
    %     res =  compute_switch_triggered_average(cellid,...
    %             'post',2,'exclude_final',0, 'condition_residual',0,...
    %           'include_str','true(size(data.trials.hit == 1))',...
    %           'which_switch',which_switch, 'n_shuffles', 250,...
    %           'save_file',1,'mask_other_switch',0,'force',1);
end
if p.plot_residual
    STR_left    = res.STR_left_real;
    STR_right   = res.STR_right_real;
else
    STR_left    = res.STA_left_real;
    STR_right   = res.STA_right_real;
end
lags        = res.lags - p.lag;
%%


ax = axes;
hold(ax,'on')

posdex = lags > -0.001 & lags < max_t;
negdex = lags <  0.001 & lags > min_t;
sat = .7; % .9
satdiff = 0; % .2
ax.TickDir = 'out';
%    clrs(2,:) = [0 .2 .8];
%    clrs(end-1,:) =[0 .8 .2];


shadedErrorBar(lags(negdex),nanmean(STR_left(:,negdex)),nansem(STR_left(:,negdex)),{'color', dp.right_color,'parent',ax},[],sat)
shadedErrorBar(lags(negdex),nanmean(STR_right(:,negdex)),nansem(STR_right(:,negdex)),{'color', dp.left_color},[],sat)
shadedErrorBar(lags(posdex),nanmean(STR_left(:,posdex)),nansem(STR_left(:,posdex)),{'color', dp.left_color},[],sat-satdiff)
shadedErrorBar(lags(posdex),nanmean(STR_right(:,posdex)),nansem(STR_right(:,posdex)),{'color', dp.right_color},[],sat-satdiff)
%%
ylabel('\Delta Firing Rate (Hz)')
xlabel(['Time from ' which_switch ' state change (s)'])
%pbaspect([2 1 1])
plot([0 0], ylim, 'k--')

ax.TickDir = 'out';

%%

sig = nan(size(res.pval));
sig(abs(res.pval-.5) > .45) = 1;
sigside = res.pval > .5;
posdex  = (lags > -0.001)';
negdex  = (lags <  0.001)';
maxy = max(ylim())*.95;
lw = 5;
pre_left = sigside==1 & negdex;
pre_right = sigside==0 & negdex;
post_left = sigside==1 & posdex;
post_right = sigside==0 & posdex;
plot(lags(pre_left),sig(pre_left).*maxy,'color',dp.left_color,'linewidth',lw)
plot(lags(pre_right),sig(pre_right).*maxy,'color',dp.right_color,'linewidth',lw)
plot(lags(post_left),sig(post_left).*maxy,'color',dp.right_color,'linewidth',lw)
plot(lags(post_right),sig(post_right).*maxy,'color',dp.left_color,'linewidth',lw)
xlim([min_t max_t])

title(['Cell ' num2str(cellid)])

if p.save_fig
    fig_name = [ num2str(res.cellid) '_' which_switch '_STA'];
    print(gcf,  fullfile(dp.fig_dir,fig_name),fig_type,'-painters')
end

