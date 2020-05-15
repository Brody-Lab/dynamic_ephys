% READ ME TYLER: This script loads the array_data, vec_data for each cell using a summary database I made. You should be able to just load array_data and vec_data using package_dyn_pbubs or whatever. And then use the included functions. Its possible get_behavior_data adds some fields to array_data that you need, but I dont think so. 
%
% I commented out, but did not delete functions that I dont think you need. 

% Set Workspace
close all; clear all;
disp('Analyzing cells!')
%cd /home/alex/Dropbox/spikes/bin

% Set Path
%addpath ~/ratter/ExperPort/MySQLUtility
%addpath ~/ratter/Analysis/Pbups
%addpath ~/ratter/Analysis/helpers
%addpath ~/ratter/ExperPort
%addpath ~/ratter/ExperPort/Analysis
%addpath ~/ratter/ExperPort/Analysis/SameDifferent/
%addpath ~/ratter/ExperPort/Utility
%addpath ~/ratter/ExperPort/bin
%addpath ~/ratter/svn_papers/TimHanks/PBupsPhys/Code/Carlosbin
%datapath = '/home/alex/Dropbox/spikes/data/';
%figpath  = '/home/alex/Dropbox/spikes/figures/';

% get list of cells
%load('../data/all_cells.mat')
%rats    = unique(allspikes.ratname);
%cells   = unique(allspikes.cellid);
%sessids = unique(allspikes.sessid);

% Setup parameters
p.reload                = 0;    % load firing data 
p.remove_initial_choice = 1;    % ignore initial model choice
p.clear_bad_strengths   = 1;    % ignore weak model state changes
p.strength_window       = 0.1;  % window for calculating model strength
p.bad_strength          = 0.05; % criteria for weak model state change
p.plot_strength         = 1;    % plot model state strength
p.plot_ignore           = 1;    % plot ignored state changes
p.firing_map_plot_trajectories = 0;
p.eval_dt               = 1e-3;

% Iterate over cells for each cell, get all the info.
for i=1:length(allspikes.cellid)
    % Load trial and spike data
%    p.ratname = allspikes.ratname{i};
%    [array_data, vec_data, cellid, sessid] = get_behavior_data(datapath, allspikes.cellid(i), allspikes.sessid(i), p);
    % GET array_data, and vec_data in some manner....

    % load just mean model trajectory for each trial
    load(['../model/model_mean_' num2str(allspikes.sessid(i)) '.mat'])
    model_mean = model_mean(vec_data.good);
    
    % compute model analysis
    array_data = compute_model_state(array_data, model_mean);
    array_data = quantify_model_state_strength(array_data,p);
    array_data = smooth_model_state(array_data,p);

    trial_browser(array_data,p,1);
    close all
end




