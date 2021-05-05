p.reload                = 0;    % load firing data 
p.remove_initial_choice = 1;    % ignore initial model choice
p.clear_bad_strengths   = 1;    % ignore weak model state changes
p.strength_window       = 0.1;  % window for calculating model strength
p.bad_strength          = 0; % criteria for weak model state change
p.plot_strength         = 1;    % plot model state strength
p.plot_ignore           = 1;    % plot ignored state changes
p.firing_map_plot_trajectories = 0;
p.eval_dt               = 1e-3;
p.model_smooth_wdw      = 100;
p.fit_line              = 1;
p.exclude_final         = 0;
p.final_only            = 0;
p.min_pre_dur           = 0;
p.min_post_dur          = 0;
p.which_switch          = 'model' ;
cellid   = 18181;
cellid = 16898;
which_trials = [6, 15, 17, 34, 65, 81];
[array_data, vec_data, this_sessid, this_rat] = package_dyn_phys(cellid);
[switch_to_0, switch_to_1, array_data, vec_data] = ...
    get_switches(cellid, ...
    'array_data',array_data,'vec_data',vec_data,...
    'which_switch',p.which_switch,...
    'clear_bad_strengths', p.clear_bad_strengths, ...
    'bad_strength', p.bad_strength, 'fit_line', p.fit_line,...
    'exclude_final', p.exclude_final, 'final_only', p.final_only,...
    'min_pre_dur',p.min_pre_dur,'min_post_dur',p.min_post_dur,...
    'model_smooth_wdw', p.model_smooth_wdw);
%%
which_trials = [10 13 15 17 29 32 80]
model_switch_browser(array_data,p,which_trials);
%%
model_switch_browser(array_data,p)