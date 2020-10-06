function build_set_ephys(rats)

if nargin < 1
    % Set which rats to analyze
%p.rats  = {'H037','H066', 'H084','H129','H140'};
%p.rats  = {'H130', 'H140','H141', 'H176','H190','H191'};
    p.rats  = {'H129'};
end
    
dp = set_dyn_path;

% Set Filenames
p.dynamic_dropbox     = dp.project_dir;
p.check_rats_dir      = fullfile(p.dynamic_dropbox, 'check_rats');
p.figpath_root        = fullfile(p.check_rats_dir, 'rat_figures/');
p.summary             = '_population';
p.population          = 'ephys';
p.modelfit            = 'model_fits';
p.model_data          = '/model_data';

% Set Rat datapath
p.datapath_root       = fullfile(p.dynamic_dropbox,'check_rats/ephys_data/');



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
p.model_color       = 'm';
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

% Do rat level analysis
for i=1:length(p.rats)
    [S] = load_data(p.rats{i},p);
    
    % Determine list of sessions to include in cosyne analysis
    analyze_rat_inclusion(p.rats{i},S,p);
    disp(' ')
end


keyboard
%% Do population level analysis
p.figpath_num         = '/rat_1_num_trials';
p.figpath_acc         = '/rat_1_accuracy';
p.figpath_vio         = '/rat_1_violations';
p.figpath_haz         = '/rat_1_hazard';
p.figpath_haz_bar     = '/rat_1_hazard_barrier';
p.figpath_clkr        = '/rat_1_click_rates';
p.figpath_durs        = '/rat_1_durations';
p.figpath_psych       = '/rat_1_psychometric';
p.figpath_psych_acc   = '/rat_1_psychometric_acc';
p.figpath_chrono_dur  = '/rat_1_chronometric_dur';
p.figpath_chrono      = '/rat_1_chronometric_gen';
p.figpath_excess      = '/rat_1_excess';
p.figpath_flips       = '/rat_1_flips';
p.figpath_manuscript  = fullfile(p.dynamic_dropbox,'/etc/manuscript/raw_figures/');

disp('Completed Individual Level Analysis, moving to population summaries')
disp(' ')

% Loop over pages
for i=1:p.numPages
    % Hazard Rate Independent Computations
    analyze_num_trials_population(p,i);
    analyze_accuracy_population(p,i);
    analyze_violations_population(p,i);
    analyze_hazard_population(p,i);
    analyze_hazard_barrier_population(p,i);
    analyze_click_rates_population(p,i);
    analyze_trial_durations_population(p,i);

    % Hazard Rate Computations
    % Population analysis methods will plot all hazard files in p.haz_compare
    analyze_psychometric_population(p,i);
    analyze_psychometric_acc_population(p,i);
    analyze_chronometric_population(p,i);
    analyze_excess_rates_population(p,i);
end

disp('Completed Population Level Analysis, moving to summary plots')
disp(' ')

% Compile list of rats to include in summary analysis
p.included = analyze_rat_inclusion_population(p);

% Compute summary analyses
summary_psychometric(p);
summary_psychometric_acc(p);
summary_chronometric(p);
summary_excess_rates(p);
summary_timescales(p);
summary_excess_rates_duration(p);

% Check for rats that aren't training
disp('Completed Summary Plots, Checking for rats off training')
disp(' ')
detect_rat_failure_population(p);

% Make Figures for paper
if 0
make_figure_chronometric(p)
make_figure_psychometric(p)
make_figure_excess_population(p)
end
