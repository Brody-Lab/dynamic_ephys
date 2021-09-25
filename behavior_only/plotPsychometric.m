function [ax, res] = plotPsychometric(ttype, outcome, varargin)
% [b edges ax] = plotPsychometric(ttype, outcome, varargin)
% ttype     ntrialsx1 regressor - e.g. click difference
% outcome   ntrialsx1 logical - e.g. animal went right

p = inputParser;
addParamValue(p,'fitLineColor',[.85 .85 .95]);
addParamValue(p,'lineColor',[0 0 .6]);
addParamValue(p,'dataColor', [0 0 0]);
addParamValue(p,'dataFaceColor', [0 0 0]);
addParamValue(p,'dataMarkerSize', 18);
addParamValue(p,'dataLineStyle', '.');
addParamValue(p,'dataLineWidth', 1);
addParamValue(p,'dataShaded',0);
addParamValue(p,'fitLineStyle', '-');
addParamValue(p,'fitLineWidth', 1);
addParamValue(p,'figHandle', []);
addParamValue(p,'axHandle', []);
addParamValue(p,'edges', []);
addParamValue(p,'nbin', 8);
addParamValue(p,'grid', 1);
addParamValue(p,'errorbar', 'binomial');
addParamValue(p,'compute_fit',true);
addParamValue(p,'plotfit', true);
addParamValue(p,'plotdata', true);
addParamValue(p,'ploterrorbar', true);
addParamValue(p,'errorbarshaded', false);
addParamValue(p,'use_fmincon', false);
addParamValue(p,'min_trials_to_plot', 10);
addParamValue(p,'fit_lapse', true);
parse(p,varargin{:});
nbin    = p.Results.nbin;
edges   = p.Results.edges;
b       = [];
if ~isempty(p.Results.axHandle)
    ax=p.Results.axHandle;
else
    fig = figure;
    ax = axes;
end

if isempty(edges) & isempty(nbin)
    ctrs        = unique(ttype);
    edges       = unique(ttype);
    [~,~,bins]  = unique(ttype);
elseif isempty(edges)
    maxmag      = max(abs(ttype));
    minttype    = min(ttype);
    maxttype    = max(ttype);
    if minttype < 0 & maxttype > 0 
        edges   = linspace(-maxmag,maxmag,nbin+1);
    else
        edges   = linspace(min(ttype),max(ttype),nbin+1);
    end
    ctrs            = (edges(1:end-1)+edges(2:end))/2;
    [bins, edges]   = discretize(ttype,edges);
else
    ctrs            = (edges(1:end-1)+edges(2:end))/2;
    [bins, edges]   = discretize(ttype,edges);
end

smoothDifficulty = min(ttype):.1:max(ttype);

% get the betas for a psychometric function
if p.Results.compute_fit
    [b, r, j, covb, mse] = fit_four_param_psycho(ttype,outcome,p.Results.use_fmincon,...
        'fit_lapse',p.Results.fit_lapse);
    % evaluate the fitted psychometric on many points
    %yfit = fourParamPsychometric(b,smoothDifficulty);
    res.b    = b(:);
    res.r    = r(:);
    res.mse  = mse(:);
    res.covb = covb;
    [yfit, delta]   = nlpredci(@fourParamPsychometric,smoothDifficulty,...
        b, r, 'covar', covb);
end

outcome = outcome(~isnan(bins));
bins    = bins(~isnan(bins));

% bin the data
nWentRight  = accumarray(bins, outcome, [], @sum);
n           = accumarray(bins, outcome, [], @length);

% compute error bars
if strcmp(p.Results.errorbar,'binomial')
    [mu, pci]   = binofit(nWentRight,n);
    up          = mu-pci(:,1);
    dwn         = pci(:,2)-mu;
elseif p.Results.ploterrorbar & strcmp(p.Results.errorbar,'gaussian')
    mu  = accumarray(bins, outcome, [], @nanmean);
    ci  = 1.96*accumarray(bins, outcome, [], @nansem);
    up  = ci;
    dwn = ci;
else
    mu  = accumarray(bins, outcome, [], @nanmean);
    up  = [];
    dwn = [];
    ci  = [];
end
hold(ax,'on')

if p.Results.compute_fit & p.Results.plotfit
    H   = shadedErrorBar(smoothDifficulty, yfit, [delta' delta'], ...
        {'color',p.Results.fitLineColor}, 1);
    
%     plot(ax,smoothDifficulty, yfit, p.Results.fitLineStyle,...
%         'color',p.Results.fitLineColor,'linewidth',p.Results.fitLineWidth);
end

if p.Results.plotdata
    mu_plot = mu;
       
    if p.Results.ploterrorbar
        up_plot     = up;
        dwn_plot    = dwn;
        nan_ind     = n < p.Results.min_trials_to_plot;
        mu_plot(nan_ind)    = nan;
        up_plot(nan_ind)    = nan;
        dwn_plot(nan_ind)   = nan;
        
        if p.Results.dataShaded
            H = shadedErrorBar(ctrs(1:length(mu)), mu_plot, [up_plot dwn_plot],...
                {'color',p.Results.dataColor}, 1);
            H.edge(1).Annotation.LegendInformation.IconDisplayStyle = 'off';
            H.edge(2).Annotation.LegendInformation.IconDisplayStyle = 'off';
        else
        plt = errorbar(ax,ctrs(1:length(mu)),mu_plot,up_plot,dwn_plot,...
            p.Results.dataLineStyle,'color',p.Results.dataColor,...
            'markerfacecolor',p.Results.dataFaceColor,...
            'linewidth',p.Results.dataLineWidth,...
            'markersize',p.Results.dataMarkerSize,'capsize',0);
        end
    else
        this_ctrs = ctrs(1:length(mu));
        plt = plot(ax,this_ctrs(:),mu_plot(:),...
            p.Results.dataLineStyle,'color',p.Results.dataColor,...
            'markerfacecolor',p.Results.dataFaceColor,...
            'linewidth',p.Results.dataLineWidth,...
            'markersize',p.Results.dataMarkerSize);
    end
end

if max(outcome) == 1
    ylim(ax,[0 1]);
end
xlim(ax, edges([1 end]));
if p.Results.grid
    ax.YGrid = 'on';
end

res.n       = n(:);
res.mu      = mu(:);
res.edges   = edges(:);
res.ctrs    = ctrs(:);
res.up      = up(:);
res.dwn     = dwn(:);