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
p.ratlist           = ratlist;
tuning_curves_dir   = fullfile(spikes_bin_dir, 'tuning_curves');
ratter_dir          = '~/ratter';
tim_code_dir        = fullfile(ratter_dir, 'Manuscripts/TimHanks/PBupsPhys/Code/');
p.code_dir          = code_dir;
p.project_dir       = project_dir;
p.check_rats_dir    = fullfile(p.project_dir, 'check_rats');
p.check_rats_figdir = fullfile(p.check_rats_dir, 'rat_figures/');
p.behav_data_dir    = fullfile(project_dir, 'data');
p.data_dir          = p.behav_data_dir;
p.example_data_dir  = fullfile(fileparts( mfilename('fullpath')),'example_data');
p.behav_data_filefn = @(rat) fullfile(behav_data_dir, [rat '.mat']);
p.behav_incl_filefn = @(rat, haz) fullfile(behav_data_dir, [rat '_inclusion' num2str(haz) '.mat']);
p.tuning_curve_dir  = tuning_curves_dir;
p.model_fits_dir    = fullfile(project_dir, 'results');
p.model_mean_dir    =  fullfile(project_dir,'model_mean');
p.model_dir         =  fullfile(project_dir,'model');
p.spikes_dir        = spikes_dir;
p.fig_dir           = fig_dir;
p.figpath_root      = fig_dir;
p.figpath_psych_acc = 'psycho_';
p.model_color       = [1 .6 .6];
p.psth_fig_dir      = fullfile(fig_dir, 'PSTH');
p.sta_fig_dir       = fullfile(fig_dir, 'STA');
p.sta_dir           = fullfile(spikes_dir, 'sta');
p.celldat_filename  = celldat_filename;
p.celldat_dir       = celldat_dir;
p.spikes_bin_dir    = spikes_bin_dir;
p.tuning_curves_dir = tuning_curves_dir;
p.ratter_dir        = ratter_dir;
p.tim_code_dir      = tim_code_dir;
p.ephys_summary_dir = fullfile(spikes_bin_dir, 'ephys_summary');
p.population        = 'ephys';
% parameters for which sessions to include
p.haz               = 1;
p.hazStr            = strrep(num2str(p.haz),'.','p');
p.haz_compare       = [1];
p.nt                = 10;
p.nt_max            = 500000;
p.nt_min            = 10;
p.ns                = 1000;
p.clearViolations   = 0;

% Parameters for Figures
p.color             = 'b';
p.compare_color     = ['b', 'r'];
p.optimal_color     = 'k';
p.model_color       = [255 140 140]/255.;
cm2                 = color_set(2);
p.left_color        = cm2(1,:);
p.right_color        = cm2(2,:);
p.pref_color  = [.8 .25 .8];
p.npref_color = [.8 .65 .25];

%p.model_color       = 'm';
p.nice_color        = {[0 140 54]./255, [48 127 255]./255 };
p.savefigure        = 1;
p.figcols           = 3;
p.figrows           = 4;
p.figwidth          = 8.5;
p.figheight         = 11;
p.figfontsize       = 12;
p.figtitlesize      = 12;
p.figaxissize       = 12;
p.singlewidth       = 4;
p.singleheight      = 4;
p.singlefontsize    = 16;
p.singletitlesize   = 16;
p.singleaxissize    = 16;
p.figvisible        = 0;
p.markersize        = 2;
p.errorbarsize      = 1;
p.mfontsize         = 20;
p.mheight           = 6;
p.mwidth            = 6;

% Parameters for Excess Click Rates
p.excess.dt         = 0.0001;
p.excess.wStd       = 500;
p.excess.ww         = p.excess.wStd*4;
p.excess.win_dur    = 0.5; 
p.excess.colorR     = 'rm';
p.excess.colorL     = 'bc';
p.excess.mean       = 1;
p.excess.plotOptimal= 0;
p.excess.normalize  = 1;
p.excess.n_by_max   = 0;
p.excess.lim        = 4*(1/p.excess.win_dur);
p.excess.plotFit    = 0;

% Parameters for Rat Inclusion
p.include.nt        = 50;
p.include.acc       = .65;
p.include.vio       = .35;%.25;
p.include.r1        = 39;
p.include.haz       = p.haz;
p.include.hazBar    = 0;
p.include.minT      = 1;
p.include.maxT      = 1.5;
p.include.ns        = 1;
p.include.nt_total  = 1000;
p.include.save      = 1;

% Parameters for optimal linear agent
p.optimal.compute   = 0;
p.optimal.lambda    = -4.1;
p.optimal.dt        = 0.0001;
p.optimal.plotNoiseOptimal = 0;

% Parameters for model behavior
p.model.compute     = 0;
p.model.plot        = 0;


p.fw = 2.5;
p.msz = 12;
p.fsz = 10.5;
set(0, 'defaultaxesfontsize',p.fsz);
set(0,'defaultaxeslinewidth',1)
set(0,'DefaultAxesTitleFontWeight','normal');


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