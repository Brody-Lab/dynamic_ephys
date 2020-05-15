
% Set Workspace
close all; clear all;
disp('Analyzing cells!')
cd /home/alex/Dropbox/spikes/bin

% Set Path
addpath ~/ratter/ExperPort/MySQLUtility
addpath ~/ratter/Analysis/Pbups
addpath ~/ratter/Analysis/helpers
addpath ~/ratter/ExperPort
addpath ~/ratter/ExperPort/Analysis
addpath ~/ratter/ExperPort/Analysis/SameDifferent/
addpath ~/ratter/ExperPort/Utility
addpath ~/ratter/ExperPort/bin
addpath ~/ratter/Manuscripts/TimHanks/PBupsPhys/Code/Carlosbin
datapath = '/home/alex/Dropbox/spikes/data/';
figpath  = '/home/alex/Dropbox/spikes/figures/';

% get list of cells
cell_list = dyn_cells_db;   % That has to be run once to create cell_list
cells = cell2mat(extracting(cell_list, 'cellid'));


% Setup parameters
p.reload                = 0;    % load firing data 
p.remove_initial_choice = 1;    % ignore initial model choice
p.clear_bad_strengths   = 1;    % ignore weak model state changes
p.strength_window       = 0.1;  % window for calculating model strength
p.bad_strength          = 0; % criteria for weak model state change
p.plot_strength         = 1;    % plot model state strength
p.plot_ignore           = 1;    % plot ignored state changes
p.firing_map_plot_trajectories = 0;
p.eval_dt               = 1e-3;
p.model_smooth_wdw      = 100;
p.fit_line              = 1;

to0 =[];
to1 = [];

% for each cell, get all the info.
for i=1:10:length(cells)
    disp(i)
    % Load trial and spike data
    [array_data, vec_data, this_sessid, this_rat] = package_dyn_phys(cells(i));
    p.ratname = this_rat;

    % throw out bad trials from array_data
    array_data = cleanup_array_data(array_data, vec_data);

    % add generative state switches to array_data
    array_data = compute_state_switches(array_data);
    
    % load just mean model trajectory for each trial
    load(['../model/model_mean_' num2str(this_sessid) '.mat'])
    model_mean = model_mean(vec_data.good);
    for j = 1:length(model_mean)
        model_mean(j).mean = movmean(model_mean(j).mean, p.model_smooth_wdw);
    end
  
    array_data = compute_model_state(array_data, model_mean);
    array_data = quantify_model_state_strength(array_data,p);
    array_data = smooth_model_state(array_data,p);
    array_data = compute_gen_state(array_data);
    to0 =[to0 [array_data.model_switch_to_0_strength]];
    to1 =[to1 [array_data.model_switch_to_1_strength]];
    trial_browser(array_data,p,1);
    close all
end

figure;hold on
histogram([-to0 to1],500)
ylabel('count')
xlabel('Slope of state change')
set(gca,'fontsize',16)
pbaspect([2 1 1])

ten = quantile([-to0 to1], .1)
plot([ten ten], ylim(), 'r--')

f = quantile([-to0 to1], .05)
plot([f f], ylim(), 'm--')


f = quantile([-to0 to1], .0225)
plot([f f], ylim(), 'c--')

save('strength_of_state_changes.mat','to1','to1')

