function [Rex, Lex, Rstd, Lstd, fit_tc,boot] = build_excess_rates(rate_one, rate_zero, choice,p,timepoint);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rdex = choice;                      % Vector with right choice trials
Ldex = ~choice;                     % Vector with left choice trials
Rrate1 = (1/p.excess.dt).*rate_one(:,Rdex);  % Divide click rates into left and right choice
Rrate0 = (1/p.excess.dt).*rate_zero(:,Rdex); % Divide click rates into left and right choice
Lrate1 = (1/p.excess.dt).*rate_one(:,Ldex);
Lrate0 = (1/p.excess.dt).*rate_zero(:,Ldex);
Rrate = [Rrate1 Rrate0];
Lrate = [Lrate1 Lrate0];
numperBinR = sqrt(sum(isnan(Rrate),2));
numperBinL = sqrt(sum(isnan(Lrate),2));
if size(Rrate,2) ~=0
    Rex  = nanmean(Rrate,2);
    Rstd = nanstd(Rrate,0,2)./numperBinR;
else
    Rex  = [];
    Rstd = [];
end
if size(Lrate,2) ~=0
    Lex  = nanmean(Lrate,2);
    Lstd = nanstd(Lrate,0,2)./numperBinL;
else
    Lex  = [];
    Lstd = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trim last data point
%%%%%% DOWNSAMPLING CHANGE
%Rex   = Rex(1:end-1);
%Lex   = Lex(1:end-1);
%Rstd  = Rstd(1:end-1);
%Lstd  = Lstd(1:end-1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalize Curves, and fit timescale
if length(Rex) > 1 &&  length(Lex) > 1
%%%%%% DOWNSAMPLING CHANGE
%    r = sum(Rex).*p.excess.dt;
%    l = sum(Rex).*p.excess.dt;
    r = sum(Rex).*p.excess.dt.*p.excess.ds;
    l = sum(Rex).*p.excess.dt.*p.excess.ds;

    Rex   = Rex./r;
    Lex   = Lex./l;
    Rstd  = Rstd./r;
    Lstd  = Lstd./l;
    dex = length(p.excess.dt:p.excess.dt*p.excess.ds:p.excess.win_dur);
    fit_tc = fit(timepoint(end-dex+1:end)', Rex(end-dex+1:end), 'exp1');
%    fit_tc = fit(timepoint(p.excess.wStd*3:end)', Rex(p.excess.wStd*3:end), 'exp1');
    bs = zeros(100,1);
    for i=1:100
        Rboot = Rex + randn(size(Rex)).*Rstd;
        fit_boot = fit(timepoint(end-dex+1:end)', Rboot(end-dex+1:end), 'exp1');
        bs(i) = fit_boot.b;
    end
    boot.b = mean(bs);
    boot.std = std(bs);
else
    fit_tc = [];    
    boot   = [];
end





