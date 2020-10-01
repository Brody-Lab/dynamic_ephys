
% Set Workspace
close all; clear all;
disp('Analyzing cells!')
dp = set_dyn_path;
datapath = dp.spikes_dir;
figpath = dp.fig_dir;
allcellspath = fullfile(datapath, 'all_cells.mat');
model_dir = dp.model_dir;

% get list of cells
load(allcellspath)
rats    = unique(allspikes.ratname);
cells   = unique(allspikes.cellid);
sessids = unique(allspikes.sessid);

% Setup parameters
p.reload                = 0;    % load firing data 
p.remove_initial_choice = 1;    % ignore initial model choice
p.clear_bad_strengths   = 0;    % ignore weak model state changes
p.strength_window       = 0.1;  % window for calculating model strength
p.bad_strength          = 0; % criteria for weak model state change
p.plot_strength         = 1;    % plot model state strength
p.plot_ignore           = 1;    % plot ignored state changes
p.firing_map_plot_trajectories = 0;
p.eval_dt               = 1e-3;
% model smooth window = 100

% for each cell, get all the info.
for i=1:length(allspikes.cellid)
    % Load trial and spike data
    p.ratname = allspikes.ratname{i};
    [array_data, vec_data, cellid, sessid] = get_behavior_data(datapath, allspikes.cellid(i), allspikes.sessid(i), p);

    % load just mean model trajectory for each trial
    load(fullfile(model_dir, ['model_mean_' num2str(allspikes.sessid(i)) '.mat']))
    model_mean = model_mean(vec_data.good);
    
    % compute model analysis
    for j = 1:length(array_data)
        model_mean(j).mean = movmean(model_mean(j).mean,50);
    end
    array_data = compute_model_state(array_data, model_mean);
    array_data = quantify_model_state_strength(array_data,p);
    array_data = smooth_model_state(array_data,p);

    trial_browser(array_data,p,1);
    close all
end




