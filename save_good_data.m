function [S] = save_good_data(rat, S, p)

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

p.behav_data_filefn(rat)
% filename = [p.datapath_root rat '/inclusion_' strrep(num2str(p.haz),'.','p')];
% save(filename, 'INCLUDE', 'session', 'trialCount','numSession','rat','p','c')
% disp(['Data saved to: ' filename]);

n_done = S.n_done(c.include);
[S] = load_session_data_(session,n_done,p);
p.population = [p.population '_' strrep(num2str(p.haz),'.','p') '_'];
if INCLUDE & p.include.save
    disp(['Saving data formatted for running Bings model'])
save_data(rat,S,p)
end




function [S] = load_session_data_(sessids,n_done,p)


% Pull data 
S = get_sessdata(sessids);
if ~isempty(sessids)
    for j=1:length(sessids)
        if sum(isnan(S.pd{j}.hits))>0
            S.pd{j}.hits(isnan(S.pd{j}.hits)) = 0;
        end
        S.acc(j) =  sum(S.pd{j}.hits) / length(S.pd{j}.hits);
        S.vio(j) =  sum(S.pd{j}.violations) / length(S.pd{j}.violations);
        S.haz(j) =  S.pd{j}.bupsdata{end}.hazard;
        S.r1(j)  =  S.pd{j}.bupsdata{end}.rate_1;
        S.r2(j)  =  S.pd{j}.bupsdata{end}.rate_2;
        S.pd{j}.stimdata = []; % Clear stim data
        S.n_done(j) = n_done(j);
    end
    S.peh = [];
end


if p.clearViolations
    for i=1:length(sessids)
        % need to document change in analysis...excess click rates changing????
        % need to separate out different types of violations?
        dex =  logical(S.pd{i}.violations);
        S.pd{i}.hits(dex)       = [];
        S.pd{i}.violations(dex) = [];
        S.pd{i}.sounds(dex)     = [];
        S.pd{i}.samples(dex)    = [];
        S.pd{i}.memory_gaps(dex)= [];
        S.pd{i}.sides(dex)      = [];
        S.pd{i}.rt_task(dex)    = [];
        S.pd{i}.stims(dex)      = [];
        S.pd{i}.bupsdata        = {S.pd{i}.bupsdata{~dex}}';
        S.pd{i}.n_left(dex)     = [];
        S.pd{i}.n_right(dex)    = [];
    end
else
    for i=1:length(sessids)
        S.pd{i}.bupsdata = {S.pd{i}.bupsdata{1:end-1}}';
    end
end