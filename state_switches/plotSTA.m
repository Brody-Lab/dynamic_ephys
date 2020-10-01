function plotSTA(cellid, which_switch, savefig)
    if nargin < 2
        which_switch = 'model';
    end
     res =  compute_switch_triggered_average(cellid,'post',2,'which_switch',which_switch, 'n_shuffles', 1000,'save_file',1,'mask_other_switch',1);
    STR_left    = res.STR_left_real;
    STR_right   = res.STR_right_real;
    lags        = res.lags;

    figure;
    hold on
    clrs = color_set(6);
    posdex = lags > -0.001;
    negdex = lags <  0.001;
    sat = .7; % .9
    satdiff = 0; % .2
%    clrs(2,:) = [0 .2 .8];
%    clrs(end-1,:) =[0 .8 .2];   
 
    shadedErrorBar(lags(negdex),nanmean(STR_left(:,negdex)),nansem(STR_left(:,negdex)),{'color', clrs(1,:)},[],sat)
    shadedErrorBar(lags(negdex),nanmean(STR_right(:,negdex)),nansem(STR_right(:,negdex)),{'color', clrs(end,:)},[],sat)
    shadedErrorBar(lags(posdex),nanmean(STR_left(:,posdex)),nansem(STR_left(:,posdex)),{'color', clrs(end-1,:)},[],sat-satdiff)
    shadedErrorBar(lags(posdex),nanmean(STR_right(:,posdex)),nansem(STR_right(:,posdex)),{'color', clrs(2,:)},[],sat-satdiff)
    xlim([-.5 1])
    set(gca, 'fontsize',16)
    ylabel('\Delta Firing Rate (Hz)')
    xlabel(['Time from ' which_switch ' state change (s)'])
    pbaspect([2 1 1])
    plot([0 0], ylim, 'k--')


sig = nan(size(res.pval));
sig(abs(res.pval-.5) > .45) = 1;
sigside = res.pval > .5;
posdex  = (res.lags > -0.001)';
negdex  = (res.lags <  0.001)';
maxy = min(ylim())*.95;
plot(res.lags(sigside==1 & negdex),sig(sigside==1 & negdex).*maxy,'color',clrs(end,:),'linewidth',10)
plot(res.lags(sigside==0 & negdex),sig(sigside==0 & negdex).*maxy,'color',clrs(1,:),'linewidth',10)
plot(res.lags(sigside==1 & posdex),sig(sigside==1 & posdex).*maxy,'color',clrs(2,:),'linewidth',10)
plot(res.lags(sigside==0 & posdex),sig(sigside==0 & posdex).*maxy,'color',clrs(end-1,:),'linewidth',10)

    title(['Cell ' num2str(cellid)])
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 6];
    fig.PaperPositionMode = 'Manual';
    fig.PaperSize = [6 6];
    if nargin > 2 && savefig
        print(gcf,  ['../figures/STA/' num2str(res.cellid) '_' which_switch '_STA'],'-dsvg')
    end

