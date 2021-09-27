function [res fh ax] = plotSTA(cellid, varargin)
p = inputParser;
addParameter(p,'which_switch', 'model')
addParameter(p,'save_fig', 1)
addParameter(p,'fig_num', [])
addParameter(p,'plot_residual',1)
addParameter(p,'lag',0)
addParameter(p,'max_t',.55)
addParameter(p,'min_t',-.55)
addParameter(p,'fig_type','svg')
addParameter(p,'ylims',[])
addParameter(p,'n_shuffles',250)
addParameter(p,'recompute',0)
addParameter(p, 'clear_bad_strengths',1 )
addParameter(p, 'bad_strength', 0)
addParameter(p, 't_buffers', [0 0])
addParameter(p, 'min_post_dur', 0)
addParameter(p, 'min_pre_dur', 0)
addParameter(p, 'alpha', 0.05)
parse(p,varargin{:})
p = p.Results;
which_switch    = p.which_switch;
save_fig        = p.save_fig;
max_t = p.max_t;
min_t = p.min_t;
fig_type = ['-d' p.fig_type];

if isempty(p.fig_num)
    fh = figure; clf;
    fpos = [4 4];
else
    fh = figure(p.fig_num); clf
    fpos = get(fh,'position');
end
set(fh,'position',[fpos([1 2]) 6 3],'papersize',[6 3],'paperpositionmode','auto')

dp = set_dyn_path;
if isstruct(cellid)
    res = cellid;
    cellid = res.cellid;
    which_switch = res.which_switch;
else
    res =  compute_switch_triggered_average(cellid, 'force', p.recompute,...
        'post',2,...
        'which_switch',which_switch, 'n_shuffles', p.n_shuffles,...
        'save_file',1,'mask_other_switch',1,...
        'bad_strength', p.bad_strength, ...
        't_buffers', p.t_buffers, ...
        'min_pre_dur', p.min_pre_dur,'min_post_dur', p.min_post_dur,...
        'clear_bad_strengths',p.clear_bad_strengths);
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
shadedErrorBar(lags(negdex),nanmean(STR_left(:,negdex)),...
    nansem(STR_left(:,negdex)),{'color', dp.right_color,'parent',ax},[],sat)
shadedErrorBar(lags(negdex),nanmean(STR_right(:,negdex)),...
    nansem(STR_right(:,negdex)),{'color', dp.left_color},[],sat)
shadedErrorBar(lags(posdex),nanmean(STR_left(:,posdex)),...
    nansem(STR_left(:,posdex)),{'color', dp.left_color},[],sat-satdiff)
shadedErrorBar(lags(posdex),nanmean(STR_right(:,posdex)),...
    nansem(STR_right(:,posdex)),{'color', dp.right_color},[],sat-satdiff)
%%
ylabel('\Delta Firing Rate (Hz)')
xlabel(['time from ' which_switch ' state change (s)'])
%pbaspect([2 1 1])
if isempty(p.ylims)
    ylims = ylim;
else
    ylims = p.ylims;
end
ylim(ylims)
plot([0 0], ylims, 'k--')

ax.TickDir = 'out';

%%
alpha = p.alpha;
sig = nan(size(res.pval));
sig(abs(res.pval-.5) > (.5-alpha/2)) = 1;
sigside = res.pval > .5;
posdex  = (lags > -0.001)';
negdex  = (lags <  0.001)';
maxy = max(ylims)*.95;
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

title(['Cell ' num2str(cellid)],'fontweight','normal')

if p.save_fig
    fig_name = [ num2str(res.cellid) '_' which_switch '_STA'];
    print(gcf,  fullfile(dp.fig_dir,fig_name),fig_type,'-painters')
end

