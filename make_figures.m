% this code regenerates the panels from the figures in Boyd-Meredith &
% Piet et al

% Before running this code, you should
% open set_dyn_path.m and edit the path configuration fields
edit set_dyn_path

% you also need to make sure to init and update the accumulation-model
% submodule. you can do this by typing `git submodule init` followed by
% `git submodule update` into the command line from the directory where
% you've cloned this repo

%% run path config
dp = set_dyn_path;
%% make fig 1 panels and S1-4 panels
print_behavior_panels;
%% make fig 2 panels and S5-6 panels
close all
clear 
print_ephys_panels
%% make fig 3 panels
close all
clear 
print_tuning_curve_panels;
%% make fig 4 panels
close all
clear 
print_state_switch_panels;
%% make fig 5 panels
close all 
clear 
use_switches = 1;
print_tuning_curve_panels;
%%
close all
clear 
%% extra supplementary stuff
% fig s1-4 get made in print_behavior_panels
% fig s5-6 get made in print_ephys_panels
%% fig s7 
random_walk_check

%% fig s8
posterior_dist_supp
%% fig s9
posterior_model_validation

%% fig s10 rank 1 review response
rank1_review_response

%% fig s11
% this was used to produce the supplementary figure on state changes
% and is useful for looking at different sets of inclusion parameters
check_switch_inclusion;
%%
% strength_histograms;