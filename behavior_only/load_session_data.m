function [S] = load_session_data(sessids)


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
        S.n_done(j) = length(S.pd{j}.hits);
    end
    S.peh = [];
end

for i=1:length(sessids)
    S.pd{i}.bupsdata = {S.pd{i}.bupsdata{1:end-1}}';
end


