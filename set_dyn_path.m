celldat_dir = '~/Dropbox/spikes/cell_packager_data';
spikes_bin_dir  = '~/Dropbox/spikes/bin';
tuning_curves_dir = '~/Dropbox/spikes/bin/tuning_curves/';
ratter_dir  = '~/ratter';

tim_code_dir = fullfile(ratter_dir, '~/ratter/Manuscripts/TimHanks/PBupsPhys/Code/');

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