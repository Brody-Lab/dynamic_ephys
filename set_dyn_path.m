function p = set_dyn_path(pathup, user)
if nargin < 1
    pathup = false;
end

if nargin < 2
    user    = 'Tyler';
end

switch user
    case 'Tyler'
        project_dir         = '~/projects/pbups_dyn/';
        spikes_dir          = fullfile(project_dir, 'data/phys');
        fig_dir             = fullfile(project_dir,'figures');
        code_dir            = fullfile(project_dir, 'code');

    case 'Alex'
        project_dir         = '~/Dropbox/spikes/';
        spikes_dir          = project_dir;
        fig_dir             = fullfile(project_dir,'figures');

    otherwise
        fprintf('you should pick a spikes directory')
        keyboard
end
celldat_filename    = 'dyn_multi_db.mat';
celldat_dir         = fullfile(spikes_dir, 'cell_packager_data');
spikes_bin_dir      = fullfile(spikes_dir, 'bin');

tuning_curves_dir   = fullfile(spikes_bin_dir, 'tuning_curves');
ratter_dir          = '~/ratter';
tim_code_dir        = fullfile(ratter_dir, 'Manuscripts/TimHanks/PBupsPhys/Code/');
p.behav_data_dir    = fullfile(project_dir, 'data');
p.tuning_curve_dir  = tuning_curves_dir;
p.model_fits_dir    = fullfile(project_dir, 'results');
p.model_mean_dir    =  fullfile(project_dir,'model_mean');
p.model_dir    =  fullfile(project_dir,'model');
p.spikes_dir        = spikes_dir;
p.fig_dir           = fig_dir;
p.psth_fig_dir      = fullfile(fig_dir, 'PSTH');
p.sta_fig_dir      = fullfile(fig_dir, 'STA');
p.celldat_filename  = celldat_filename;
p.celldat_dir       = celldat_dir;
p.spikes_bin_dir    = spikes_bin_dir;
p.tuning_curves_dir = tuning_curves_dir;
p.ratter_dir        = ratter_dir;
p.tim_code_dir      = tim_code_dir;
p.ephys_summary_dir = fullfile(spikes_bin_dir, 'ephys_summary');


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
    
    addpath(genpath(code_dir));

end


if ~exist(p.psth_fig_dir)
    mkdir(p.psth_fig_dir);
end