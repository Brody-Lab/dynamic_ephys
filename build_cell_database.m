
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
addpath ~/ratter/svn_papers/TimHanks/PBupsPhys/Code/Carlosbin
datapath = '/home/alex/Dropbox/spikes/data/';
figpath = '/home/alex/Dropbox/spikes/figures/';

% Setup parameters
p.reload    = 1;
p.normalize = 1;
p.raster    = 1;
p.sort_by_state = -1;
p.click_triggered = 1;
p.sort_by_state_dur = 'none';
p.firing_map_plot_trajectories = 0;
 
rats = {'H037', 'H066', 'H084', 'H129', 'H140'};
allspikes = struct();
for i=1:length(rats)
    p.ratname   = rats{i};
    disp(p.ratname)
    % get spike data
    spikeinfo   = get_spikeinfo(p.ratname);
    spikeinfo.ratname = cell(length(spikeinfo.sessid),1);
    [spikeinfo.ratname{:}] = deal(p.ratname);
    if i==1
        allspikes = spikeinfo;
    else
        allspikes   = merge_spikeinfo(allspikes, spikeinfo);
    end
end

% save spike info
save([datapath 'all_cells.mat'], 'allspikes')

% for each cell, get all the info.
for i=1:length(allspikes.cellid)
    p.ratname = allspikes.ratname{i};
    [array_data, vec_data, cellid, sessid] = get_behavior_data(datapath, allspikes.cellid(i), allspikes.sessid(i), p);

    % compute average firing rate and statistics during trials
    [mean_fr, var_fr, base_fr, dprime, last_mean_fr, last_var_fr, last_dprime]  = compute_average_firing_rate(array_data,vec_data);
    allspikes.mean_fr(i)    = mean_fr;
    allspikes.var_fr(i)     = var_fr;
    allspikes.dprime(i)     = dprime;
    allspikes.fano(i)       = mean_fr/var_fr;
    allspikes.base_fr(i)    = base_fr;
    allspikes.last_mean_fr(i)   = last_mean_fr;
    allspikes.last_var_fr(i)    = last_var_fr;
    allspikes.last_dprime(i)    = last_dprime;
    allspikes.last_fano(i)      = last_mean_fr/last_var_fr;   

    % compute time-dependent statistics
    % compute average firing rate split on choice
    % compute dprime as function of time
    % compute fano-factor as function of time
    % get time-indexes for which each cell becomes side selective. bootroc.m
end

% resave database with statistics
save([datapath 'all_cells.mat'], 'allspikes')

dex = allspikes.mean_fr > 5;
% 40% cells have firing rate > 5 hz
% 15% of cells have dprime > 0.1
% 9% have both
figure; 
subplot(2,3,1); hold on;
plot(abs(allspikes.dprime(dex)), abs(allspikes.last_dprime(dex)), 'bo')
plot([0 0.6],[0 0.6], 'k--')
ylabel('200ms d''')
xlabel('full trial d''')
axis square 

subplot(2,3,2); hold on;
plot(allspikes.fano(dex), allspikes.last_fano(dex), 'bo')
plot([0 2],[0 2], 'k--')
ylabel('200ms fano factor''')
xlabel('full trial fano factor''')
axis square 

subplot(2,3,3); hold on;
plot(allspikes.mean_fr(dex), allspikes.last_mean_fr(dex), 'bo')
plot([0 100],[0 100], 'k--')
ylim([0 100])
xlim([0 100])
axis square 
ylabel('200ms f.r.')
xlabel('full trial f.r.')

subplot(2,3,4); hold on;
plot(allspikes.var_fr(dex), allspikes.last_var_fr(dex), 'bo')
plot([0 100],[0 100], 'k--')
ylim([0 100])
xlim([0 100])
axis square 
ylabel('200ms variance')
xlabel('full trial variance')

subplot(2,3,5); hold on;
plot(allspikes.base_fr(dex), allspikes.mean_fr(dex), 'bo')
plot([0 100],[0 100], 'k--')
axis square 
ylabel('trial f.r.')
xlabel('pre trial f.r.')
ylim([0 100])
xlim([0 100])









