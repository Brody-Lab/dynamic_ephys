function [] = plot_sta_dprime(res, fh)
if nargin < 2
    fh = figure;
end
clrs = color_set(6);

sig = nan(size(res.pval));
sig(abs(res.pval-.5) > .45) = 1;
%supersig = nan(size(res.pval));
%supersig(abs(res.pval-.5) > .49) = 1;
sigside = res.pval > .5;
posdex  = (res.lags > -0.001)';
negdex  = (res.lags <  0.001)';

plot(res.lags, res.dprime_shuff, 'color', [1 1 1].*.7)

hold on
%%plot(res.lags, res.pval,'c'); %%a% debug
plot(res.lags, res.dprime_real, 'color', [.65 .05 .05], 'linewidth', 2);
xlabel(['Time from ' res.params.which_switch ' state change (s)'])
ylabel('d''')
set(gca, 'fontsize',16)
plot(xlim, [0 0], 'k')
plot( [0 0],ylim, 'k')
xlim([-.5 1])
maxy = 1.1.*max(abs(res.dprime_real(res.lags > -.5 & res.lags < 1)));
ylim([-1 1].*maxy)

plot(res.lags(sigside==1 & negdex),sig(sigside==1 & negdex).*maxy,'s','markeredgecolor',clrs(end,:),'markerfacecolor',clrs(end,:),'markersize',5)
plot(res.lags(sigside==0 & negdex),sig(sigside==0 & negdex).*maxy,'s','markeredgecolor',clrs(1,:), 'markerfacecolor',clrs(1,:),'markersize',5)
plot(res.lags(sigside==1 & posdex),sig(sigside==1 & posdex).*maxy,'s','markeredgecolor',clrs(2,:), 'markerfacecolor',clrs(2,:),'markersize',5)
plot(res.lags(sigside==0 & posdex),sig(sigside==0 & posdex).*maxy,'s','markeredgecolor',clrs(end-1,:), 'markerfacecolor',clrs(end-1,:),'markersize',5)

plot(res.lags(sigside==1 & negdex),sig(sigside==1 & negdex).*maxy,'color',clrs(end,:),'linewidth',5)
plot(res.lags(sigside==0 & negdex),sig(sigside==0 & negdex).*maxy,'color',clrs(1,:),'linewidth',5)
plot(res.lags(sigside==1 & posdex),sig(sigside==1 & posdex).*maxy,'color',clrs(2,:),'linewidth',5)
plot(res.lags(sigside==0 & posdex),sig(sigside==0 & posdex).*maxy,'color',clrs(end-1,:),'linewidth',5)
box off
set(gcf, 'color', 'w')
title(['Cell ' num2str(res.cellid)])
pbaspect([2 1 1])


