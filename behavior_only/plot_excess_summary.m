cd ~/Dropbox/spikes/bin/behavior_only/
load ~/Dropbox/dynamic/model_fits/ANALYSIS/ephys/group_analysis.mat

for i=1:length(F)
    close all
    if isfield(F{i}, 'clicks')
        plot_excess_rates(F{i}.clicks)
        print(['/home/alex/Dropbox/spikes/figures/task_figure/' F{i}.rat '_excess'], '-dsvg')
    else
        addpath ~/Dropbox/dynamic/check_rats
        load(['~/Dropbox/dynamic/model_fits/FITS/ephys/fit_analytical_' F{i}.rat '.mat'])
        p = struct();
        clicks = compute_excess_rates_data(data,p) 
    end
end


