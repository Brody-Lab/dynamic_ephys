function p = set_dyn_path(pathup, user)
if nargin < 1
    pathup = false;
end

if nargin < 2
    user    = 'Tyler';
end

switch user
    case 'Tyler'
        spikes_dir          = '~/projects/pbups_dyn/data/phys';
    case 'Alex'
        spikes_dir          = '~/Dropbox/spikes';
    otherwise
        fprintf('you should pick a spikes directory')
        keyboard
end
spikes_fig_dir      = '~/Dropbox/spikes/figures';

celldat_filename    = 'dyn_multi_db.mat';
celldat_dir         = fullfile(spikes_dir, 'cell_packager_data');
spikes_bin_dir      = fullfile(spikes_dir, 'bin');
tuning_curves_dir   = fullfile(spikes_bin_dir, 'tuning_curves');
ratter_dir          = '~/ratter';
tim_code_dir        = fullfile(ratter_dir, 'Manuscripts/TimHanks/PBupsPhys/Code/');
p.psth_fig_dir      = fullfile(spikes_fig_dir, 'PSTH');
p.celldat_filename  = celldat_filename;
p.celldat_dir       = celldat_dir;
p.spikes_bin_dir    = spikes_bin_dir;
p.tuning_curves_dir = tuning_curves_dir;
p.ratter_dir        = ratter_dir;
p.tim_code_dir      = tim_code_dir;

if pathup
    addpath(tim_code_dir, fullfile(tim_code_dir, 'Carlosbin'))
    addpath(fullfile(ratter_dir, 'ExperPort/bin'))
    addpath(fullfile(ratter_dir, 'ExperPort/MySQLUtility'))
    
    addpath(fullfile(ratter_dir, 'ExperPort/Analysis'))
    addpath(fullfile(ratter_dir, 'ExperPort/SameDifferent'))
    addpath(fullfile(ratter_dir, 'ExperPort/HandleParam'))
    
    addpath(fullfile(ratter_dir, 'Analysis/Pbups'))
    addpath(fullfile(ratter_dir, 'Analysis/helpers'))
    
    addpath(celldat_dir)
    addpath(spikes_bin_dir)
    addpath(tuning_curves_dir)
end


