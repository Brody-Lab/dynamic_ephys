function [fig, ax] = plot_excess_rates(clicks,varargin)

p = inputParser;
addParameter(p, 'plot_model', 0);
addParameter(p,'fig_num',[]);
addParameter(p,'left_color',[]);
addParameter(p,'right_color',[]);
parse(p,varargin{:});
plot_model = p.Results.plot_model;

if ~isempty(p.Results.fig_num)
    fig = figure(p.Results.fig_num);
else
    fig = figure;
end
ax = axes;
smoothing = 25;
%nice_color        = {[0 140 54]./255, [48 127 255]./255 };

hold on;
%plot(clicks.timepoint, clicks.RexRat, 'color', p.nice_color{1});
%plot(clicks.timepoint, clicks.LexRat, 'color', p.nice_color{2});
R = movmean(clicks.RexRat,smoothing);
L = movmean(clicks.LexRat,smoothing);
Rs= movmean(clicks.RstdRat,smoothing);
Ls= movmean(clicks.LstdRat,smoothing);

shadedErrorBar(clicks.timepoint,R,Rs, {'Color', p.Results.right_color},0)  
shadedErrorBar(clicks.timepoint,L,Ls, {'Color', p.Results.left_color},0)  

if plot_model & isfield(clicks, 'RexN')
Rn = movmean(clicks.RexN,smoothing);
Ln = movmean(clicks.LexN,smoothing);
Rsn= movmean(clicks.RstdN,smoothing);
Lsn= movmean(clicks.LstdN,smoothing);


shadedErrorBar(clicks.timepoint,Rn,Rsn,'m',0)  
shadedErrorBar(clicks.timepoint,Ln,Lsn, 'm',0)  
end

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



