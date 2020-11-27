function [clicks] = compute_excess_rates_data(data,p)

    p.excess.win_dur    = 0.5;
    p.excess.wStd       = 500;%80;
    p.excess.ds         = 10;
    p.excess.dt         = 0.0001;
    p.excess.ww         = p.excess.wStd*4;
    p.excess.mean       = 1;

    p.excess.normalize  = 1;
    p.excess.minT       = p.excess.win_dur + .05;
    p.excess.window=exp(-.5.*(1:p.excess.ww).^2./p.excess.wStd^2)./sum(exp(-.5.*(1:p.excess.ww).^2./p.excess.wStd^2));
    p.excess.full_steps = length(p.excess.dt:p.excess.dt:p.excess.minT);
    p.excess.steps      = length(p.excess.dt:p.excess.dt*p.excess.ds:p.excess.minT);
    p.compute_model = 0;
    p.optimal.compute = 0;
    p.optimal.lambda = -4.1;
    p.nt_min = 1;

    [rate_one, rate_zero, nt, choice, sides,dex] = build_rates_data(data,p);




% Set up output variable
c.p         = p.excess;
c.noise     = p.optimal;
%%% DOWNSAMPLING CODE CHANGE
%c.timepoint = -p.excess.minT + p.excess.dt:p.excess.dt:-p.excess.dt;
c.timepoint = -p.excess.minT + p.excess.dt:p.excess.dt*p.excess.ds:-p.excess.dt;
c.xlim      = [-p.excess.minT+p.excess.dt+p.excess.wStd*p.excess.dt -p.excess.dt];
c.nt        = nt;

%regress, normalize, and fit timeconstant
t = c.timepoint;
[c.RexRat, c.LexRat, c.RstdRat, c.LstdRat, c.fitRat,c.bootRat] = build_excess_rates(rate_one, rate_zero, choice,p,t);
[c.RexOpt, c.LexOpt, c.RstdOpt, c.LstdOpt, c.fitOpt,c.bootOpt] = build_excess_rates(rate_one, rate_zero, sides, p,t);

clicks = c;
