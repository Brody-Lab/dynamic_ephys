function [] = plot_parameter_comparisons(F,group,dyn_boxplot);



% Plotting boundaries
bounds = [-12 4;-.1 30; -.1 10; -.1 30; -.1 50;  -.1 1; -.1 .7; -1 1; -.1 0.5];

% load Bing parameters
load ~/Dropbox/model_fits/FITS/best_fits_v35_rats.mat;
B = allfits;
B(:,3) = B(:,3)./40;
Bse = allci;
Bse(:,3) = Bse(:,3)./40;

% parameter names
param_names = { '\lambda',  '\sigma_a^2', '\sigma_s^2', '\sigma_i^2','B', '\phi', '\tau_{\phi}', 'bias', 'lapse'};

% Remove parameters
if strcmp(group, 'with_sticky_bounds')
    error('not implemented')

elseif strcmp(group, 'a_x_stream') | strcmp(group, 'a_fits_3') | ...
        strcmp(group, 'dataset_2') | strcmp(group, 'BING') | ...
        strcmp(group, 'BING_noL') | strcmp(group, 'dataset_2_noL') | ...
        strcmp(group, 'dataset_2_prior')| ...
        strcmp(group, 'dataset_2_gaussian_prior') | ...
        strcmp(group, 'BING_ALICE') | strcmp(group, 'ephys')
    B(:,5) = [];   
    param_names(5) = [];
    bounds(5,:) = [];

elseif strcmp(group, 'a_no_si')
    B(:,5) = [];   
    B(:,4) = [];   
    Bse(:,4) = [];
    param_names(5) = [];
    param_names(4) = [];
    bounds(5,:) = [];
    bounds(4,:) = [];

elseif strcmp(group, 'a_no_sa')
    B(:,5) = [];   
    B(:,4) = [];   
    B(:,2) = [];   
    
    Bse(:,4) = [];
    Bse(:,2) = [];
    param_names(5) = [];
    param_names(4) = [];
    param_names(2) = [];
    bounds(5,:) = [];
    bounds(4,:) = [];
    bounds(2,:) = [];   
else
error('What model verison is this? Check group name')
end

H = [];
P = [];
AR = [];
AB = [];


%%F(8) = [];

for i=1:length(param_names)
    [h,p,avg_rat,avg_bing] = plot_parameter_dist(F,i,...
        param_names(i),B(:,i),Bse(:,i),bounds(i,:),group,...
        'fig_num',i,'point_plot',0);
    H=[H h];
    P=[P p];
    AR=[AR avg_rat];
    AB=[AB avg_bing];
end


