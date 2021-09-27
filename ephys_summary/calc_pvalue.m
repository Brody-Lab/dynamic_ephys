function p = calc_pvalue(frate, conds, varargin)
% function p = calc_pvalue(frate, conds, varargin)
%
% Calculates p-value of difference in firing rate between conditions at
% each time step.
%
% INPUT:
%
%   frate:  trial x time firing rate matrix
%   conds:  condition for each trial
%
% OUTPUT:
% 
%   p:      p-value as a function of time
%

% stats = 'anova';        % statistical method to compare different conditions





%% Default Parameters
p = inputParser;
addParameter(p, 'stats', 'regression')
parse(p,varargin{:});
stats = p.Results.stats;

cond_list = unique(conds);
% perform statistical test at each time point
p = nan(1,size(frate,2));
switch stats
    case 'anova'
        for i=1:size(frate,2)
            if sum(~isnan(frate(:,i)))>1
                p(i) = anova1(frate(:,i), conds, 'off');
            else
                p(i) = nan;
            end
        end
    case 'ttest'
        [h, p, ci] = ttest2(frate(conds==cond_list(1),:), frate(conds==cond_list(2),:));
    case 'roc'
        [auc, p] = slidingROC(frate(conds==cond_list(1),:), frate(conds==cond_list(2),:));
    case 'regression'
        for i=1:size(frate,2)
            [b, bint, ~, ~, rstats] = regress(frate(:,i), [ones(size(conds)) conds]);
            p(i) = rstats(3);
        end
    otherwise
        error('invalid stats test in calc_pvalue');
end