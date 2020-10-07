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

