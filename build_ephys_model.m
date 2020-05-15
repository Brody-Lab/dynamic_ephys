close all;
clear all;
% a script for getting the backwards pass model output for ephys rats
cd ~/Dropbox/spikes/bin/

addpath ~/Dropbox/dynamic/model_fits/analytical_model/
addpath ~/ratter/ExperPort/MySQLUtility
addpath ~/ratter/Analysis/Pbups
addpath ~/ratter/Analysis/helpers
addpath ~/ratter/ExperPort
addpath ~/ratter/ExperPort/Analysis
addpath ~/ratter/ExperPort/Analysis/SameDifferent/
addpath ~/ratter/ExperPort/Utility
addpath ~/ratter/ExperPort/bin
addpath ~/ratter/Manuscripts/TimHanks/PBupsPhys/Code/Carlosbin

% get database of all cells
load('../data/all_cells.mat')

% list of rats
rats = unique(allspikes.ratname);

for i=1:length(rats)
disp(rats{i})    
% get index of which sessions to analyze
sessids = unique(allspikes.sessid(ismember(allspikes.ratname, rats{i})));

for j=1:length(sessids)
    disp(sessids(j))

    % load data
    S = load_session_data(sessids(j));
    data = format_data(S);

    % load params
    load(['/home/alex/Dropbox/dynamic/model_fits/ANALYSIS/ephys/fit_analysis_analytical' rats{i} '.mat'])
    params = fit.final;

    % compute model output for each trial
    p_in.error_tolerance    = 1e-4;
    p_in.compute_dist       = 1;
    p_in.da_grid            = 0.1;
    p_in.da                 = 0.1;
    p_in.avals              = -10:p_in.da:10;
    p_in.just_pdf           = 1;
    [model,model_p] = accumulation_model(data,params,'posterior',p_in);

    % save
    model_p.sessid = sessids(j);
    model_p.ratname = rats{i};
    s = whos('model');
    if s.bytes > 2e9
        save(['../model/model_posterior_' num2str(sessids(j)) '.mat'], '-v7.3', 'model','model_p')
    else
        save(['../model/model_posterior_' num2str(sessids(j)) '.mat'], 'model','model_p')
    end

    % save separate file with just mean because its more compact and sufficient for some analyses
    for k=1:length(model)
        model_mean(k).mean  = model(k).posterior.mean;
        model_mean(k).T     = model(k).posterior.T;
    end
    save(['../model/model_mean_' num2str(sessids(j)) '.mat'], 'model_mean','model_p')
end
end


