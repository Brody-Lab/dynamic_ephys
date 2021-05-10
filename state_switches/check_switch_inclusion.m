p.reload                = 0;    % load firing data 
p.remove_initial_choice = 1;    % ignore initial model choice
p.clear_bad_strengths   = 1;    % ignore weak model state changes
p.strength_window       = 0.1;  % window for calculating model strength
p.bad_strength          = 0; % criteria for weak model state change
p.plot_strength         = 1;    % plot model state strength
p.plot_ignore           = 1;    % plot ignored state changes
p.firing_map_plot_trajectories = 0;
p.eval_dt               = 1e-3;
p.model_smooth_wdw      = 150;
p.fit_line              = 1;
p.exclude_final         = 0;
p.final_only            = 0;
p.min_pre_dur           = 0;
p.min_post_dur          = 0;

p.which_switch          = 'model' ;
cellid   = 18181;
cellid = 16898;


p.remove_initial_choice = 0; %1;
p.change_bounds = [-1 1].*0;

[array_data_, vec_data_, this_sessid, this_rat] = package_dyn_phys(cellid);

[switch_to_0, switch_to_1, array_data, vec_data] = ...
    get_switches(cellid, ...
    'array_data',array_data_,'vec_data',vec_data_,...
    'which_switch',p.which_switch,...
    'clear_bad_strengths', p.clear_bad_strengths, ...
    'bad_strength', p.bad_strength, 'fit_line', p.fit_line,...
    'exclude_final', p.exclude_final, 'final_only', p.final_only,...
    'min_pre_dur',p.min_pre_dur,'min_post_dur',p.min_post_dur,...
    'model_smooth_wdw', p.model_smooth_wdw,...
    'change_bounds',p.change_bounds,...
    'remove_initial_choice',p.remove_initial_choice);
%%
smooth_level = [1 100 250 500];
figure(1); clf
s = [];
dur_x = 0:.1:2;
rate_x = 0:.2:5;

for jj = 1:length(smooth_level)
    %%
    [switch_to_0, switch_to_1, array_data, vec_data] = ...
        get_switches(cellid, ...
        'array_data',array_data_,'vec_data',vec_data_,...
        'which_switch',p.which_switch,...
        'clear_bad_strengths', p.clear_bad_strengths, ...
        'bad_strength', p.bad_strength, 'fit_line', p.fit_line,...
        'exclude_final', p.exclude_final, 'final_only', p.final_only,...
        'min_pre_dur',p.min_pre_dur,'min_post_dur',p.min_post_dur,...
        'model_smooth_wdw', smooth_level(jj),...
        'change_bounds',p.change_bounds,...
        'remove_initial_choice',p.remove_initial_choice);
    
    m_state_durs = [array_data(:).model_state_durs];
    m_switch_rate = [array_data(:).model_switches_per_time];
    g_state_durs = [array_data(:).gen_state_durs];
    g_switch_rate = [array_data(:).gen_switches_per_time];
    
    s(1,jj) = subplot(2,length(smooth_level),jj);
    h = histogram(m_state_durs,'binedges',dur_x)
    hold on
    histogram(g_state_durs,'displaystyle','stairs','binedges',dur_x)

    title(sprintf('moving average %i ms',smooth_level(jj)))
    xlabel('duration (s)')
    box off
    
    
    s(2,jj) = subplot(2,length(smooth_level),length(smooth_level)+jj);
    plot([ 1 1],[0 200],'k')
    hold on
    h=histogram(m_switch_rate,'binedges',rate_x)
    hold on
    histogram(g_switch_rate,'displaystyle','stairs','binedges',rate_x)

    box off
    xlabel('switches / s')  
    
    drawnow

end
linkaxes(s(1,:))
linkaxes(s(2,:))
%%
disp(1)
p.remove_initial_choice = 0 ;
p.change_bounds = [-1 1].*0;
p.bad_strength = 0;
p.model_smooth_wdw = 100;
p.t_buffers = [.15 .15];
p.min_pre_dur = 0;
p.min_post_dur = 0;
[switch_to_0, switch_to_1, array_data, vec_data] = ...
    get_switches(cellid, ...
    'array_data',array_data_,'vec_data',vec_data_,...
    'which_switch',p.which_switch,...
    'clear_bad_strengths', p.clear_bad_strengths, ...
    'bad_strength', p.bad_strength, 'fit_line', p.fit_line,...
    'exclude_final', p.exclude_final, 'final_only', p.final_only,...
    'min_pre_dur',p.min_pre_dur,'min_post_dur',p.min_post_dur,...
    'model_smooth_wdw', p.model_smooth_wdw,...
    'change_bounds',p.change_bounds,...
    'remove_initial_choice',p.remove_initial_choice,...
    't_buffers',p.t_buffers);
%%

bad_trials = [9 15 95]
model_switch_browser(array_data,p,bad_trials);
%%
which_trials = [10 13 15 17 29 32 80]
which_trials = [6, 15, 17, 34, 65, 81];
which_trials = [6, 15, 17, 34, 65, 81 10 13 15 17 29 32 80]
model_switch_browser(array_data,p,which_trials);
%%
model_switch_browser(array_data,p)
%%
