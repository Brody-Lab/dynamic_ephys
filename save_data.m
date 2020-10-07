function [] = save_data(rat, S, p);

save_name_include = fullfile(p.behav_data_dir, [rat '.mat']);
save_name_full = fullfile(p.behav_data_dir, [rat '_full_dataset.mat']);

% Format for Bing's model
data = format_data(S);
disp(['#trials: ' num2str(length(data))]);
% Save
if ~exist(p.behav_data_dir,'dir')
    mkdir(p.behav_data_dir)
end
save(save_name_include, 'data', 'p', 'rat');
disp(['Data saved to:  ' save_name_include]);

save(save_name_full, 'S', 'p', 'rat');
disp(['Data saved to: ' save_name_full]);


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