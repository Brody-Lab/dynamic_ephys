function res = plotSTA(cellid, which_switch, savefig)
    if nargin < 2
        which_switch = 'model';
    end
    dp = set_dyn_path;
%         res =  compute_switch_triggered_average(cellid,'post',2,...
%           'which_switch',which_switch, 'n_shuffles', 1000,...
%           'save_file',1,'mask_other_switch',1);
    res =  compute_switch_triggered_average(cellid,...
            'post',2,'exclude_final',0, 'condition_residual',0,...
          'include_str','true(size(data.trials.hit == 1))',...
          'which_switch',which_switch, 'n_shuffles', 250,...
          'save_file',1,'mask_other_switch',0,'force',1);
    STR_left    = res.STR_left_real;
    STR_right   = res.STR_right_real;
    lags        = res.lags;
%%
    fh = figure(2); clf;
    set(fh,'position',[2 2 6 3],'papersize',[6 3],'paperpositionmode','auto')
    ax = axes;
    hold(ax,'on')
    max_t = .75;
    min_t = -.45;
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
posdex  = (res.lags > -0.001)';
negdex  = (res.lags <  0.001)';
maxy = max(ylim())*.95;
lw = 5;
plot(res.lags(sigside==1 & negdex),sig(sigside==1 & negdex).*maxy,'color',dp.left_color,'linewidth',lw)
plot(res.lags(sigside==0 & negdex),sig(sigside==0 & negdex).*maxy,'color',dp.right_color,'linewidth',lw)
plot(res.lags(sigside==1 & posdex),sig(sigside==1 & posdex).*maxy,'color',dp.right_color,'linewidth',lw)
plot(res.lags(sigside==0 & posdex),sig(sigside==0 & posdex).*maxy,'color',dp.left_color,'linewidth',lw)
xlim([min_t max_t])

    title(['Cell ' num2str(cellid)])

    if nargin > 2 && savefig
        fig_name = [ num2str(res.cellid) '_' which_switch '_STA'];
        print(gcf,  fullfile(dp.fig_dir,fig_name),'-dsvg','-painters')
    end

