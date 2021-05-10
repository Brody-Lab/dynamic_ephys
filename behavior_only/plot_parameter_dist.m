function [h,p,avg_rat, avg_bing,fh] = plot_parameter_dist(F,i,param_name,...
    B,Bse, bounds, group, varargin);
p = inputParser;
addParameter(p, 'bing_boxplot', 0);
addParameter(p, 'dyn_boxplot', 0);
addParameter(p, 'point_plot', 0);
addParameter(p, 'fig_num', 1);
parse(p,varargin{:});

dp = set_dyn_path;
bing_boxplot    = p.Results.bing_boxplot;
dyn_boxplot     = p.Results.dyn_boxplot
point_plot      = p.Results.point_plot;
fig_num         = p.Results.fig_num;

bing_color = [.8  .9 1];

bing_color = [.95 .8 .9];
bing_color = [.9 .8 .9];
dyn_color = dp.model_color; % [1 .55 .55]
if nargin < 8
    dyn_boxplot = 0;
end

plot_ratnames = 1;
no_sort = 0;
fh = figure(fig_num);  clf
plot([.9 3.6],[0 0],'-','color',[1 1 1].*.7)
box off
hold on;
set(fh,'position',[5 5 dp.fw dp.fw])


P = [];
if no_sort
    h   = (1/length(F))*.9;
    sp  = (1/length(F))*.05;
    start = (1/length(F));
    bs  = start*(length(F)+1);
    d   = .7;
    NUM = length(F);
    
    for j=1:length(F)
        if isstruct(F{j})
            lb = F{j}.final(i) - F{j}.se(i);
            width = F{j}.se(i)*2;
            if isnan(lb)
                lb = eps;
                width = eps*2;
            end
            p = F{j}.final(i);
            P = [P; p];
            if ~dyn_boxplot & ~point_plot
                rectangle('Position', [lb start*j+sp width h],...
                    'EdgeColor', [d d d],'FaceColor',[d d d]);
                plot([p p], [start*j+sp start*j+h+sp], 'k','linewidth',4)
            end
        end
    end
else
    LB      = [];
    WIDTH   = [];
    P       = [];
    ratnames = {};
    for j=1:length(F)
        if isstruct(F{j})
            lb = F{j}.final(i) - F{j}.se(i);
            width = F{j}.se(i)*2;
            p = F{j}.final(i);
            if isnan(lb)
                lb = eps;
                width = eps*2;
            end
            LB = [LB lb];
            WIDTH = [WIDTH width];
            P = [P p];
            ratnames{j} = F{j}.rat;
        end
    end
    % Sort that shit
    [P,dex] = sort(P);
    ratnames = ratnames(dex);
    LB = LB(dex);
    WIDTH = WIDTH(dex);
    WIDTH(WIDTH ==0) = eps;
    h   = (1/length(P))*.9;
    sp  = (1/length(P))*.05;
    start = (1/length(P));
    bs  = start*(length(P)+1);
    d   = .7;
    NUM = length(P);
    rat_color = [1 .8 .8];
    if ~dyn_boxplot  & ~point_plot
        rectangle('Position',[bounds(1) start*(1) (bounds(2)-bounds(1)) 1],'EdgeColor',rat_color,'FaceColor',rat_color);
        
        for j=1:length(P)
            rectangle('Position', [LB(j) start*j+sp WIDTH(j) h],'EdgeColor', [d d d],'FaceColor',[d d d]);
            plot([P(j) P(j)], [start*j+sp start*j+h+sp], 'k','linewidth',2)
            if plot_ratnames
                text(LB(j) + WIDTH(j) + mean(WIDTH), start*j+sp+h/2, ratnames{j})
            end
        end
    end
end
if dyn_boxplot
    %%
    p = median(P);
    %lb = p - 1.96*sem(B);
    lims = percentile(P, [25 75]);
    lb = lims(1);
    width = diff(lims);
    %box_y = start*(NUM+1)+Bsp+Bstart*(j-1);
    box_y = start*(j)+sp ;
    
    if 0
        plot([min(P) max(P)], [1 1].*box_y+h/2, 'k','linewidth',2)
        hold on
        rectangle('Position', [lb box_y width h],...
            'EdgeColor', [1 1 1].*0,'FaceColor','w');
        plot([p p], [box_y box_y+h], 'k','linewidth',2)
        plot([0 0], [start-sp box_y+h],'k--')
        ylim([start-sp box_y+h])
    else
        h = 1;
        box_y = 1;
        plot([1 1].*box_y+h/2, [min(P) max(P)], 'k','linewidth',2)
        hold on
        rectangle('Position', [box_y lb  h width],...
            'EdgeColor', [1 1 1].*0,'FaceColor','w');
        plot([box_y box_y+h], [p p], 'k','linewidth',2)
        plot([start-sp box_y+h], [0 0], 'k--')
        xlim([start-sp box_y+h])
    end
    
    %text(lims(2) + width/3, box_y+h/1.25, 'Dynamic rats (n=5)')
end
if point_plot
    %%
    box_y = 1;
    h = 1;
    Px = box_y+h/2+(rand(size(P))-.5)*h/2'
    if ~dyn_boxplot
        p = median(P);
        plot([box_y box_y+h], [p p], 'k','linewidth',1.5)
    end
    
    plot(Px, P, 'o',...
        'markerfacecolor',dyn_color,'markersize',9,...
        'markeredgecolor','w')
    ex_rat = 'H037';
    ex =  find(ismember(ratnames,ex_rat));
    text(Px(ex)+.1,P(ex)-.4,ex_rat)
    plot(Px(ex), P(ex), 'o',...
        'markerfacecolor',dyn_color,'markersize',9,...
        'markeredgecolor','k')
end


% Plot Bing's rats
[B,dex] = sort(B);
Bse     = Bse(dex);

% new plotting values
Bh      = (1/length(B))*.8;
Bsp     = (1/length(B))*.1;
Bstart  = (1/length(B));
Bbs     = start*(length(P)+1);

% set bing background color


if ~bing_boxplot & ~point_plot
    
    rectangle('Position',[bounds(1) start*(NUM+1) (bounds(2)-bounds(1)) 1],'EdgeColor',bing_color,'FaceColor',bing_color);
    
    % Iterate Bing rats
    for j=1:length(B)
        c       = j+NUM;
        lb      = B(j) - Bse(j);
        width   = Bse(j)*2;
        p       = B(j);
        rectangle('Position', [lb start*(NUM+1)+Bsp+Bstart*(j-1) width Bh],'EdgeColor', [d d d],'FaceColor',[d d d]);
        plot([p p], [start*(NUM+1)+Bsp+Bstart*(j-1) start*(NUM+1)+Bh+Bsp+Bstart*(j-1)], 'k','linewidth',2)
    end
    plot([0 0], [start-sp start*(NUM+1)+1],'k:')
    ylim([start-sp start*(NUM+1)+1])
    
elseif bing_boxplot
    %%
    p = median(B);
    %lb = p - 1.96*sem(B);
    lims = percentile(B, [25 75]);
    lb = lims(1);
    width = diff(lims);
    %box_y = start*(NUM+1)+Bsp+Bstart*(j-1);
    box_y = start*(j+1)+sp ;
    if 0
        plot([min(B) max(B)], [1 1].*box_y+h/2, 'k','linewidth',2)
        hold on
        rectangle('Position', [lb box_y width h],...
            'EdgeColor', [1 1 1].*0,'FaceColor','w');
        plot([p p], [box_y box_y+h], 'k','linewidth',2)
        plot([0 0], [start-sp box_y+h],'k--')
        ylim([start-sp box_y+h])
    else
        box_y = 2.5;
        h = 1;
        plot([1 1].*box_y+h/2, [min(B) max(B)], 'k','linewidth',2)
        hold on
        rectangle('Position', [box_y lb h width],...
            'EdgeColor', [1 1 1].*0,'FaceColor','w');
        plot([box_y box_y+h], [p p], 'k','linewidth',2)
        plot([start-sp box_y+h], [0 0], 'k--')
        xlim([start-sp box_y+h])
    end
    
    %text(lims(2) + width/3, box_y+h/1.25, 'Brunton rats (n=19)')
    
    %boxplot(B,'orientation','horizontal')
    %%
end
if point_plot
    %%
    box_y = 2.5;
    h = 1;
    if ~bing_boxplot
        p = median(B);
        plot([box_y box_y+h], [p p], 'k','linewidth',1.5)
    end
    plot(box_y+h/2+(rand(size(B))-.5)*h/2, B,'o',...
        'markerfacecolor',bing_color,...
        'markeredgecolor','w','markersize',9)
end
%%
%NUM = NUM +length(B);
if ~point_plot
    ylabel('Rats','fontsize',dp.fsz);
    xlabel(param_name,'fontsize',dp.fsz);
    set(gca,'fontsize',dp.fsz);
    version -release
    test = ans;
    
    if strcmp(test, '2013b')
        set(gca, 'YTick', [])
    else
        yticks([])
    end
    xlim(bounds)
    
else
    xlabel('environment','fontsize',dp.fsz);
    ylabel(param_name,'fontsize',dp.fsz);
    set(gca,'fontsize',dp.fsz);
    set(gca,'xtick',[1.5 3],'xticklabel',{'dynamic', 'stationary'})
    xtickangle(.15)
    ylim(bounds)
    xlim([.9 3.6])
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5


axis square
% fig = gcf;
% fig.PaperUnits = 'inches';
% fig.PaperPosition = [0 0 6 6];
% fig.PaperPositionMode = 'Manual';
% fig.PaperSize = [6 6];
%%%%%%%%%%%%%%%%%%%%%55
dp = set_dyn_path;
fname = fullfile(dp.fig_dir, [group '_param_comparison_' num2str(i)]);

print(fh,fname,'-dsvg')

%
% do population significance test
[h,p]= ttest2(B,P);
avg_rat = mean(P);
avg_bing = mean(B);

