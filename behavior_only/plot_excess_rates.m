function [] = plot_excess_rates(clicks,varargin)

smoothing = 25;
p.nice_color        = {[0 140 54]./255, [48 127 255]./255 };
if length(varargin) > 0
    figure(varargin{1}); 
    if length(varargin) > 1
        subplot(varargin{2}, varargin{3}, varargin{4}); 
    end
else
    figure; 
end
hold on;
%plot(clicks.timepoint, clicks.RexRat, 'color', p.nice_color{1});
%plot(clicks.timepoint, clicks.LexRat, 'color', p.nice_color{2});
R = movmean(clicks.RexRat,smoothing);
L = movmean(clicks.LexRat,smoothing);
Rs= movmean(clicks.RstdRat,smoothing);
Ls= movmean(clicks.LstdRat,smoothing);

shadedErrorBar(clicks.timepoint,R,Rs, {'Color', p.nice_color{1}},0)  
shadedErrorBar(clicks.timepoint,L,Ls, {'Color', p.nice_color{2}},0)  

if isfield(clicks, 'RexN') & false
Rn = movmean(clicks.RexN,smoothing);
Ln = movmean(clicks.LexN,smoothing);
Rsn= movmean(clicks.RstdN,smoothing);
Lsn= movmean(clicks.LstdN,smoothing);


shadedErrorBar(clicks.timepoint,Rn,Rsn,'m',0)  
shadedErrorBar(clicks.timepoint,Ln,Lsn, 'm',0)  
end

%plot(clicks.timepoint, clicks.fitRat.a.*exp(clicks.fitRat.b.*clicks.timepoint),'r');

xlim(clicks.xlim)
set(gca,'fontsize',20);
ylabel('stimulus weighting','fontsize',20)
xlabel('time from end of trial (s)', 'fontsize',20)
axis square
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 6];
fig.PaperPositionMode = 'Manual';
fig.PaperSize = [6 6];
ylim([-8 8])



