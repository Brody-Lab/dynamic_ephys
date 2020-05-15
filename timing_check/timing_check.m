close all
clear all

cd ~/Dropbox/spikes/bin/timing_check/

% path stuff
addpath ~/ratter/Manuscripts/TimHanks/PBupsPhys/Code/
addpath ~/ratter/Manuscripts/TimHanks/PBupsPhys/Code/Carlosbin
addpath ~/ratter/ExperPort/bin
addpath ~/ratter/Analysis/Pbups
addpath ~/ratter/ExperPort/MySQLUtility
addpath ~/ratter/ExperPort/Analysis
addpath ~/ratter/ExperPort/Analysis/SameDifferent/
addpath ~/ratter/ExperPort/HandleParam
addpath ~/ratter/Analysis/helpers
addpath ~/Dropbox/spikes/cell_packager_data
addpath ~/Dropbox/spikes/bin
addpath ~/Dropbox/spikes/bin/tuning_curves/

% load list of cells
cell_list = dyn_cells_db;   % That has to be run once to create cell_list

% For each cell, load array_data, and perform check
cellids = cell2mat(extracting(cell_list, 'cellid'));
sessids = cell2mat(extracting(cell_list, 'sessid'));
ratnames = cell2mat(extracting(cell_list, 'ratname'));

stim_aligned     = [];
not_stim_aligned = [];
aligned_names = [];
not_aligned_names = [];

for i=1:length(cellids)
    disp(i)
    cellid = cellids(i);
    [array_data vec_data sessid ratname] = package_dyn_phys(cellid);
    if vec_data.cpoke_start(1) == vec_data.stim_start(1)
        stim_aligned = [cellid; stim_aligned];
        aligned_names = unique([ratnames(i,:); aligned_names],'rows')
    else
        disp('not aligned')
        not_stim_aligned = [cellid; not_stim_aligned];
        not_aligned_names = unique([ratnames(i,:); not_aligned_names],'rows')
    end
end

save('alignments.mat','stim_aligned','not_stim_aligned','aligned_names','not_aligned_names', 'cellids','sessids');







