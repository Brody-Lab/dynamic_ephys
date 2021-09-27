function dp = set_dyn_path(pathup, user)
% function dp = set_dyn_path(pathup, user)
% example use dp = set_dyn_path(1, 'Tyler')
%
% This sets up parameters of the dynamic ephys project. File paths, figure
% colors, ratnames, inclusion criteria for the paper are all set here.
%
% A new user trying to replicate the results of the paper shoud add
% themself as a user, change the default user to their own name, and define
% the variables project_dir, spikes_dir, fig_dir, code_dir


if nargin < 1
    pathup = false;
end
if nargin < 2
    user    = 'Tyler';
end

switch user
    case 'Tyler'
        project_dir         = '~/projects/pbups_dyn/';
        code_dir            = fullfile(project_dir, 'code');
        project_dir         = '~/projects/test_pbups_dyn/';
        spikes_dir          = fullfile(project_dir, 'data/phys');
        fig_dir             = fullfile(project_dir,'figures');

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
dp.ratlist           = ratlist;
tuning_curves_dir   = fullfile(spikes_bin_dir, 'tuning_curves');
ratter_dir          = '~/ratter';
tim_code_dir        = fullfile(ratter_dir, 'Manuscripts/TimHanks/PBupsPhys/Code/');
dp.code_dir          = code_dir;
dp.project_dir       = project_dir;
dp.check_rats_dir    = fullfile(dp.project_dir, 'check_rats');
dp.check_rats_figdir = fullfile(dp.check_rats_dir, 'rat_figures/');
dp.behav_data_dir    = fullfile(project_dir, 'data');
dp.data_dir          = dp.behav_data_dir;
dp.example_data_dir  = fullfile(fileparts( mfilename('fullpath')),'example_data');
dp.behav_data_filefn = @(rat) fullfile(dp.behav_data_dir, [rat '.mat']);
dp.behav_incl_filefn = @(rat, haz) fullfile(behav_data_dir, [rat '_inclusion' num2str(haz) '.mat']);
dp.tuning_curve_dir  = tuning_curves_dir;
dp.model_fits_dir    = fullfile(project_dir, 'results');
dp.model_mean_dir    =  fullfile(project_dir,'model_mean');
dp.model_dir         =  fullfile(project_dir,'model');
dp.spikes_dir        = spikes_dir;
dp.fig_dir           = fig_dir;
dp.figpath_root      = fig_dir;
dp.figpath_psych_acc = 'psycho_';
dp.model_color       = [1 .6 .6];
dp.psth_fig_dir      = fullfile(fig_dir, 'PSTH');
dp.sta_fig_dir       = fullfile(fig_dir, 'STA');
dp.sta_dir           = fullfile(spikes_dir, 'sta');
dp.celldat_filename  = celldat_filename;
dp.celldat_dir       = celldat_dir;
dp.spikes_bin_dir    = spikes_bin_dir;
dp.tuning_curves_dir = tuning_curves_dir;
dp.ratter_dir        = ratter_dir;
dp.tim_code_dir      = tim_code_dir;
dp.ephys_summary_dir = fullfile(spikes_bin_dir, 'ephys_summary');
dp.population        = 'ephys';
% parameters for which sessions to include
dp.haz               = 1;
dp.hazStr            = strrep(num2str(dp.haz),'.','p');
dp.haz_compare       = [1];
dp.nt                = 10;
dp.nt_max            = 500000;
dp.nt_min            = 10;
dp.ns                = 1000;
dp.clearViolations   = 0;

% Parameters for Figures
dp.color             = 'b';
dp.compare_color     = ['b', 'r'];
dp.optimal_color     = 'k';
dp.model_color       = [255 140 140]/255.;
cm2                 = color_set(2);
dp.left_color        = cm2(1,:);
dp.right_color        = cm2(2,:);
dp.pref_color  = [.8 .25 .8];
dp.npref_color = [.8 .65 .25];

%p.model_color       = 'm';
dp.nice_color        = {[0 140 54]./255, [48 127 255]./255 };
dp.savefigure        = 1;
dp.figcols           = 3;
dp.figrows           = 4;
dp.figwidth          = 8.5;
dp.figheight         = 11;
dp.figfontsize       = 12;
dp.figtitlesize      = 12;
dp.figaxissize       = 12;
dp.singlewidth       = 4;
dp.singleheight      = 4;
dp.singlefontsize    = 16;
dp.singletitlesize   = 16;
dp.singleaxissize    = 16;
dp.figvisible        = 0;
dp.markersize        = 2;
dp.errorbarsize      = 1;
dp.mfontsize         = 20;
dp.mheight           = 6;
dp.mwidth            = 6;

% Parameters for Excess Click Rates
dp.excess.dt         = 0.0001;
dp.excess.wStd       = 500;
dp.excess.ww         = dp.excess.wStd*4;
dp.excess.win_dur    = 0.5; 
dp.excess.colorR     = 'rm';
dp.excess.colorL     = 'bc';
dp.excess.mean       = 1;
dp.excess.plotOptimal= 0;
dp.excess.normalize  = 1;
dp.excess.n_by_max   = 0;
dp.excess.lim        = 4*(1/dp.excess.win_dur);
dp.excess.plotFit    = 0;

% Parameters for Rat Inclusion
dp.include.nt        = 50;
dp.include.acc       = .65;
dp.include.vio       = .35;%.25;
dp.include.r1        = 39;
dp.include.haz       = dp.haz;
dp.include.hazBar    = 0;
dp.include.minT      = 1;
dp.include.maxT      = 1.5;
dp.include.ns        = 1;
dp.include.nt_total  = 1000;
dp.include.save      = 1;

% Parameters for optimal linear agent
dp.optimal.compute   = 0;
dp.optimal.lambda    = -4.1;
dp.optimal.dt        = 0.0001;
dp.optimal.plotNoiseOptimal = 0;

% Parameters for model behavior
dp.model.compute     = 0;
dp.model.plot        = 0;

dp.fw = 2.5;
dp.msz = 12;
dp.fsz = 10.5;
set(0, 'defaultaxesfontsize',dp.fsz);
set(0,'defaultaxeslinewidth',1)
set(0,'DefaultAxesTitleFontWeight','normal');

if pathup
    if 0
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
    end  
    addpath(genpath(code_dir));
end

if ~exist(dp.psth_fig_dir)
    mkdir(dp.psth_fig_dir);
end