function [fig, ax] = plot_excess_rates(clicks,varargin)

p = inputParser;
addParameter(p, 'plot_model', 0);
addParameter(p,'fig_num',[]);
addParameter(p,'left_color',[48 127 255]./255);
addParameter(p,'right_color',[0 140 54]./255);
addParameter(p,'model_color',[]);
parse(p,varargin{:});
plot_model = p.Results.plot_model;

if ~isempty(p.Results.fig_num)
    fig = figure(p.Results.fig_num);
else
    fig = figure;
end
ax = axes;
smoothing = 25;
%nice_color        = {,  };

hold on;
%plot(clicks.timepoint, clicks.RexRat, 'color', p.nice_color{1});
%plot(clicks.timepoint, clicks.LexRat, 'color', p.nice_color{2});
R = movmean(clicks.RexRat,smoothing);
L = movmean(clicks.LexRat,smoothing);
Rs= movmean(clicks.RstdRat,smoothing);
Ls= movmean(clicks.LstdRat,smoothing);

if plot_model & isfield(clicks, 'RexN')
Rn = movmean(clicks.RexN,smoothing);
Ln = movmean(clicks.LexN,smoothing);
Rsn= movmean(clicks.RstdN,smoothing);
Lsn= movmean(clicks.LstdN,smoothing);


shadedErrorBar(clicks.timepoint,Rn,Rsn,{'color',p.Results.model_color},0)  
shadedErrorBar(clicks.timepoint,Ln,Lsn, {'color',p.Results.model_color},0)  
%plot(clicks.timepoint,Rn,'color',p.Results.model_color,'linewidth',2)
%hold on
%plot(clicks.timepoint,Ln,'color',p.Results.model_color,'linewidth',2)
end

h1 = shadedErrorBar(clicks.timepoint,R,Rs, {'Color', p.Results.right_color},1)  
h2 = shadedErrorBar(clicks.timepoint,L,Ls, {'Color', p.Results.left_color},1)  
delete h1.mainLine;
delete h2.mainLine;




%plot(clicks.timepoint, clicks.fitRat.a.*exp(clicks.fitRat.b.*clicks.timepoint),'r');

xlim(clicks.xlim)
% set(gca,'fontsize',20);
% ylabel('relative weight','fontsize',20)
% xlabel('time from end of trial (s)', 'fontsize',20)
axis square
% fig = gcf;
% fig.PaperUnits = 'inches';
% fig.PaperPosition = [0 0 6 6];
% fig.PaperPositionMode = 'Manual';
% fig.PaperSize = [6 6];
ylim([-8 8])



