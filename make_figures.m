% this code regenerates the panels from the figures in Boyd-Meredith &
% Piet et al

% Before running this code, you should open set_dyn_path.m and edit the
% path configuration fields
edit set_dyn_path
%%
print_behavior_panels;
%%
close all
clear 
print_ephys_panels
%%
close all
clear 
print_tuning_curve_panels;
%%
close all
clear 
print_state_switch_panels;
%%
close all 
clear 
use_switches = 1;
print_tuning_curve_panels;
%%
% generate fig S10, which answers a reviewer's question about rank 1
% decomposition
rank1_review_response