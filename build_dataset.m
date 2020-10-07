function build_dataset(rats)

if nargin < 1
    % Set which rats to analyze
    %rats  = {'H037','H066', 'H084','H129','H140'};
    %rats  = {'H130', 'H140','H141', 'H176','H190','H191'};
    rats  = {'H129'};
end
    
p       = set_dyn_path;
p.rats  = rats;

% Do rat level analysis
for ii = 1 : length(p.rats)
    [S] = load_data(p.rats{ii},p);
    save_good_data(p.rats{ii},S,p);
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
p.figpath_manuscript  = fullfile(dp.project_dir,'/etc/manuscript/raw_figures/');

disp('Completed Individual Level Analysis, moving to population summaries')
disp(' ')

% Loop over pages
for ii=1:p.numPages
    % Hazard Rate Independent Computations
    analyze_num_trials_population(p,ii);
    analyze_accuracy_population(p,ii);
    analyze_violations_population(p,ii);
    analyze_hazard_population(p,ii);
    analyze_hazard_barrier_population(p,ii);
    analyze_click_rates_population(p,ii);
    analyze_trial_durations_population(p,ii);

    % Hazard Rate Computations
    % Population analysis methods will plot all hazard files in p.haz_compare
    analyze_psychometric_population(p,ii);
    analyze_psychometric_acc_population(p,ii);
    analyze_chronometric_population(p,ii);
    analyze_excess_rates_population(p,ii);
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
