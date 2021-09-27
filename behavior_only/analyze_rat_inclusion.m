function [] = analyze_rat_inclusion(rat, S, p);

disp(['Determining inclusion for ' rat]);

% convert gamma to largest click rate
%R1 = (40*exp(p.include.gamma))/(1+exp(p.include.gamma));

% Individual inclusion criteria
c.vio     = S.vio <= p.include.vio;
c.acc     = S.acc >= p.include.acc;
c.nt      = S.n_done >= p.include.nt;
c.haz     = S.haz == p.include.haz;
c.rates   = S.r1 <= p.include.r1;
%c.rates   = ismembertol(S.r1, R1, 1e-3);

for i=1:length(S.pd)
    if ~isempty(S.pd{i}.bupsdata)
        temp        = [S.pd{i}.bupsdata{:}];
        hbar        = temp(1).hazardBarrier;
        c.hbar(i)   = hbar <= p.include.hazBar;
        T           = [temp.T];
        minT(i)     = min(T);
        maxT(i)     = max(T);
        c.T(i)      = min(T) <= p.include.minT & max(T) >= p.include.maxT;

    else
        c.hbar(i) = 0;
        c.T(i) = 0;
    end
end


% Conjuction of all criteria
c.include   = c.vio & c.acc & c.nt & c.haz & c.rates & c.hbar & c.T;
numSession  = sum(c.include);
session     = S.sessid(c.include);
trialCount  = sum(S.n_done(c.include));

% Does rat qualify?
INCLUDE = numSession >= p.include.ns & trialCount >= p.include.nt_total;

if INCLUDE
    disp(['Rat ' rat ' pass inclusion criteria'])
else
    disp(['RAT ' rat ' FAILS INCLUSION CRITERIA'])
end

% Muscimol Infusion sessions for H061. 
if find(session == [515515, 517387, 518622])
    error('Infusion sessions!!!!')
end

% save inclusion file with session list, and inclusion yes/no
if 7~=exist([p.datapath_root rat])
    mkdir([p.datapath_root rat])
end
filename = [p.datapath_root rat '/inclusion_' strrep(num2str(p.haz),'.','p')];
save(filename, 'INCLUDE', 'session', 'trialCount','numSession','rat','p','c')
disp(['Data saved to: ' filename]);


if INCLUDE & p.include.save
    n_done = S.n_done(c.include);
    disp(['Saving data formatted for running Bings model'])
    [S] = load_session_data(session,n_done,p);
    p.population = [p.population '_' strrep(num2str(p.haz),'.','p') '_'];
    save_data(rat,S,p)
end
