function [clicks] = compute_excess_rates(S,p);

% Compute Excess Click Rates -------------------
%warning('off','all')

% Check if we have enough trials
if p.nt_min > sum(S.n_done(S.haz == p.haz))
    clicks = NaN;
    warning('on','all')
    return
end

% Build Rate Matrices, and choice vectors
if ~p.optimal.compute & ~p.compute_model
    [rate_one, rate_zero, nt, choice, sides] = build_rates(S,p);
else
    [rate_one, rate_zero, nt, choice, sides, noise_sides, model_sides] = build_rates(S,p);
end

% Set up output variable
c.p         = p.excess;
c.haz       = p.haz;
c.noise     = p.optimal;
%%% DOWNSAMPLING CODE CHANGE
%c.timepoint = -p.excess.minT + p.excess.dt:p.excess.dt:-p.excess.dt;
c.timepoint = -p.excess.minT + p.excess.dt:p.excess.dt*p.excess.ds:-p.excess.dt;
c.xlim      = [-p.excess.minT+p.excess.dt+p.excess.wStd*p.excess.dt -p.excess.dt];
c.nt        = nt;

% Regress, normalize, and fit timeconstant
t = c.timepoint;
[c.RexRat, c.LexRat, c.RstdRat, c.LstdRat, c.fitRat,c.bootRat] = build_excess_rates(rate_one, rate_zero, choice,p,t);
[c.RexOpt, c.LexOpt, c.RstdOpt, c.LstdOpt, c.fitOpt,c.bootOpt] = build_excess_rates(rate_one, rate_zero, sides, p,t);
if p.optimal.compute
    [c.RexN, c.LexN, c.RstdN, c.LstdN, c.fitN,c.bootN] = build_excess_rates(rate_one, rate_zero, noise_sides, p,t);
end
if p.compute_model
    [c.RexM, c.LexM, c.RstdM, c.LstdM, c.fitM,c.bootM] = build_excess_rates(rate_one, rate_zero, model_sides, p,t);
end

clicks = c;

%warning('on','all')




