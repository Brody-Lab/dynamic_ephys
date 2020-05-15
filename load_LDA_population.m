function [res, cellids,computed, not_computed] = load_LDA_population(which_switch, this_force)
%which_switch = 'generative'; % 'model' or 'generative'

cd ~/Dropbox/spikes/bin
addpath ~/ratter/svn_papers/TimHanks/PBupsPhys/Code/
addpath ~/ratter/svn_papers/TimHanks/PBupsPhys/Code/Carlosbin
addpath ~/ratter/ExperPort/bin
addpath ~/ratter/Analysis/Pbups
addpath ~/ratter/ExperPort/MySQLUtility
addpath ~/ratter/ExperPort/Analysis
addpath ~/ratter/ExperPort/Analysis/SameDifferent/
addpath ~/ratter/ExperPort/HandleParam
addpath ~/ratter/Analysis/helpers
addpath ~/Dropbox/spikes/cell_packager_data/
addpath ~/Dropbox/spikes/bin/tuning_curves

if nargin < 2
    this_force = 0;
end


if ~this_force && exist(['../cell_packager_data/population_LDA_' which_switch '.mat']) == 2
    load(['../cell_packager_data/population_LDA_' which_switch '.mat']);
else
    % get the list of cells
    cell_list = dyn_cells_db;
    % only look at cells with reasonable firing rates during the trial
    select_str = 'normmean > .5';
    cellids = cell2mat(extracting(cell_list, 'cellid',select_str));
    % iterate over cells and compute STAs
    computed        = zeros(size(cellids));
    not_computed    = ones(size(cellids));
    res = cell(1,length(cellids));
    nn = 1;
    cellids = [18181 16857];
     for cc = nn:length(cellids)
         try
             res{nn} =  compute_switch_triggered_average(cellids(cc),'post',2,'which_switch',which_switch, 'n_shuffles', 1000,'save_file',1,'mask_other_switch',1);

             res{nn}.STR_right_shuff = [];
             res{nn}.STR_left_shuff  = [];
             res{nn}.dprime_real     = [];
             res{nn}.dprime_shuff    = [];
             res{nn}.params          = [];
             res{nn}.pval            = [];
             disp(cc)
             nn = nn+1;
            computed(cc) = 1;
            not_computed(cc) = 0;
         catch ME
            disp(ME.message)
         end 
     end
     
    save(['../cell_packager_data/population_LDA_' which_switch '.mat'],'cellids','res','not_computed','computed')
end


