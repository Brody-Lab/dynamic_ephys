function res = model_switch_browser(array_data, p,which_trials)
dp = set_dyn_path;
if nargin < 3
    which_trials = 1:length(array_data);
end
if which_trials > length(array_data)
    disp('bad start_trial')
    which_trials = 1;
end
res = [];
figure(1); clf
hold on;
for ii=1:length(which_trials)
    
    jj = which_trials(ii);
    model_state = 2*array_data(jj).model_state-1;
    gen_state   = 2*array_data(jj).gen_state-1;
    model_left  = nan(size(model_state));
    model_right  = nan(size(model_state));
    model_left(model_state==-1) = -1;
    model_right(model_state==1) = 1;
    clf;
    hold on;
    plot(array_data(jj).model_T([1 end]),[0 0], 'color',[1 1 1].*.5)
    if isfield(array_data, 'raw_model_mean')
        plot(array_data(jj).model_T, array_data(jj).raw_model_mean, '-','color',[1 1 1].*.5);
    end
    
    
    plot(array_data(jj).model_T, array_data(jj).model_mean, '-',...
        'color',dp.model_color,'linewidth',2);
    plot(array_data(jj).model_T, model_state, ':', 'color', [.5 .5 .5])
    plot(array_data(jj).model_T, model_left, 'color', dp.left_color)
    plot(array_data(jj).model_T, model_right, 'color', dp.right_color)
    plot(array_data(jj).model_switch_to_0, ...
        zeros(1,length(array_data(jj).model_switch_to_0)), 'ro','markerfacecolor','r')
    plot(array_data(jj).model_switch_to_1, ...
        zeros(1,length(array_data(jj).model_switch_to_1)), 'ko','markerfacecolor','k')
    if isfield(array_data(jj), 'ignore_model_switch_to_0')
        plot(array_data(jj).ignore_model_switch_to_0, ...
            zeros(1,length(array_data(jj).ignore_model_switch_to_0)), 'ro')
    end
    if isfield(array_data(jj), 'ignore_model_switch_to_1')
        plot(array_data(jj).ignore_model_switch_to_1, ...
            zeros(1,length(array_data(jj).ignore_model_switch_to_1)), 'ko')
    end
    ymax = max(ylim);
    
    tvec = 1e-4:1e-4:array_data(jj).stim_end;
    gen_left = nan(size(tvec));
    gen_right = nan(size(tvec));
    gen_right(gen_state==1) = ymax+.5;
    gen_left(gen_state==-1) = ymax;
    plot(tvec, ymax+.5*(gen_state==1), ':','color',[1 1 1].*.5)
    plot(tvec, gen_right, '-','color',dp.right_color)
    plot(tvec, gen_left, '-','color',dp.left_color)
    x = ceil(max(abs(array_data(jj).model_mean))+.1) + 1;
    ylim([-x x])
    xlim([0 2])
    pbaspect([3 1 1])
    ylabel('Accumulation Value')
    xlabel('Time')
    set(gca,'fontsize',12)
    
    slopes = [];
    switch_to_0_slope = [];
    switch_to_1_slope = [];
    ignore_switch_to_0_slope = [];
    ignore_switch_to_1_slope = [];
    if p.fit_line
        
        
        if isfield(array_data(jj), 'model_switch_to_0_strength') && p.plot_strength;
            for kk=1:length(array_data(jj).model_switch_to_0)
                T = array_data(jj).model_switch_to_0(kk);
                slope = array_data(jj).model_switch_to_0_strength(kk);
                plot([T-p.strength_window, T+p.strength_window],...
                    [0-slope*p.strength_window, 0+slope*p.strength_window], 'r-')
                switch_to_0_slope(kk) = slope;
            end
            for kk=1:length(array_data(jj).model_switch_to_1)
                T = array_data(jj).model_switch_to_1(kk);
                slope = array_data(jj).model_switch_to_1_strength(kk);
                plot([T-p.strength_window, T+p.strength_window],...
                    [0-slope*p.strength_window, 0+slope*p.strength_window], 'k-')
                switch_to_1_slope(kk) = slope;
            end 
            
        end
        if isfield(array_data(jj), 'ignore_0_strength') && p.plot_strength && p.plot_ignore
            for kk=1:length(array_data(jj).ignore_model_switch_to_0)
                T = array_data(jj).ignore_model_switch_to_0(kk);
                slope = array_data(jj).ignore_0_strength(kk);
                plot([T-p.strength_window, T+p.strength_window], [0-slope*p.strength_window, 0+slope*p.strength_window], 'r--')
                ignore_switch_to_0_slope(kk) = slope;
            end
            for kk=1:length(array_data(jj).ignore_model_switch_to_1)
                T = array_data(jj).ignore_model_switch_to_1(kk);
                slope = array_data(jj).ignore_1_strength(kk);
                plot([T-p.strength_window, T+p.strength_window], [0-slope*p.strength_window, 0+slope*p.strength_window], 'k--')
                ignore_switch_to_1_slope(kk) = slope;
            end 
        end
    else
        if isfield(array_data(jj), 'model_switch_to_0_strength') ...
                && p.plot_strength;
            keyboard
            for kk=1:length(array_data(jj).model_switch_to_0)
                T = array_data(jj).model_switch_to_0(kk);
                slope = array_data(jj).model_switch_to_0_strength(kk);
                plot([T-p.strength_window, T+p.strength_window], [0-slope, 0+slope], 'r-')
            end
            for kk=1:length(array_data(jj).model_switch_to_1)
                T = array_data(jj).model_switch_to_1(kk);
                slope = array_data(jj).model_switch_to_1_strength(kk);
                plot([T-p.strength_window, T+p.strength_window], [0-slope, 0+slope], 'k-')
            end 
        end
        if isfield(array_data(jj), 'ignore_0_strength') ...
                && p.plot_strength && p.plot_ignore
            keyboard
            for kk=1:length(array_data(jj).ignore_model_switch_to_0)
                T = array_data(jj).ignore_model_switch_to_0(kk);
                slope = array_data(jj).ignore_0_strength(kk);
                plot([T-p.strength_window, T+p.strength_window], [0-slope, 0+slope], 'r--')
            end
            for kk=1:length(array_data(jj).ignore_model_switch_to_1)
                T = array_data(jj).ignore_model_switch_to_1(kk);
                slope = array_data(jj).ignore_1_strength(kk);
                plot([T-p.strength_window, T+p.strength_window], [0-slope, 0+slope], 'k--')
            end 
        end
    end
    title(['Trial ' num2str(jj)])

    if p.wait
        fprintf('press any key to continue to the next test trial')
        pause()
    else
        pause(3);
    end
    
    res(ii).T = array_data(jj).model_T;
    res(ii).gen_state = gen_state;
    res(ii).model_mean =  array_data(jj).model_mean;
    res(ii).model_state =  array_data(jj).model_state;
    res(ii).model_switch_to_0 =  array_data(jj).model_switch_to_0;
    res(ii).ignore_model_switch_to_0 =  array_data(jj).ignore_model_switch_to_0;
    res(ii).ignore_model_switch_to_1 =  array_data(jj).ignore_model_switch_to_1;
    res(ii).switch_to_0_slope = switch_to_0_slope;
    res(ii).switch_to_1_slope = switch_to_1_slope;
    res(ii).ignore_switch_to_0_slope = ignore_switch_to_0_slope;
    res(ii).ignore_switch_to_1_slope = ignore_switch_to_1_slope;
end

