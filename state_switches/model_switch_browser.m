function [] = model_switch_browser(array_data, p,which_trials)
dp = set_dyn_path;
if nargin < 3
    which_trials = 1:length(array_data);
end
if which_trials > length(array_data)
    disp('bad start_trial')
    which_trials = 1;
end

figure(1); clf
hold on;
for i=which_trials
    model_state = 2*array_data(i).model_state-1;
    gen_state   = 2*array_data(i).gen_state-1;
    model_left  = nan(size(model_state));
    model_right  = nan(size(model_state));
    model_left(model_state==-1) = -1;
    model_right(model_state==1) = 1;
    clf;
    hold on;
    plot(array_data(i).model_T([1 end]),[0 0], 'color',[1 1 1].*.5)
    plot(array_data(i).model_T, array_data(i).raw_model_mean, '-','color',[1 1 1].*.5);
    plot(array_data(i).model_T, array_data(i).model_mean, '-',...
        'color',dp.model_color,'linewidth',2);
    plot(array_data(i).model_T, model_state, ':', 'color', [.5 .5 .5])
    plot(array_data(i).model_T, model_left, 'color', dp.left_color)
    plot(array_data(i).model_T, model_right, 'color', dp.right_color)
    plot(array_data(i).model_switch_to_0, ...
        zeros(1,length(array_data(i).model_switch_to_0)), 'ro','markerfacecolor','r')
    plot(array_data(i).model_switch_to_1, ...
        zeros(1,length(array_data(i).model_switch_to_1)), 'ko','markerfacecolor','k')
    if isfield(array_data(i), 'ignore_model_switch_to_0')
        plot(array_data(i).ignore_model_switch_to_0, ...
            zeros(1,length(array_data(i).ignore_model_switch_to_0)), 'ro')
    end
    if isfield(array_data(i), 'ignore_model_switch_to_1')
        plot(array_data(i).ignore_model_switch_to_1, ...
            zeros(1,length(array_data(i).ignore_model_switch_to_1)), 'ko')
    end
    ymax = max(ylim);
    
    tvec = 1e-4:1e-4:array_data(i).stim_end;
    gen_left = nan(size(tvec));
    gen_right = nan(size(tvec));
    gen_right(gen_state==1) = ymax+.5;
    gen_left(gen_state==-1) = ymax;
    plot(tvec, ymax+.5*(gen_state==1), ':','color',[1 1 1].*.5)
    plot(tvec, gen_right, '-','color',dp.right_color)
    plot(tvec, gen_left, '-','color',dp.left_color)
    x = ceil(max(abs(array_data(i).model_mean))+.1) + 1;
    ylim([-x x])
    xlim([0 2])
    pbaspect([3 1 1])
    ylabel('Accumulation Value')
    xlabel('Time')
    set(gca,'fontsize',12)

    if p.fit_line
        if isfield(array_data(i), 'model_switch_to_0_strength') && p.plot_strength;
            for j=1:length(array_data(i).model_switch_to_0)
                T = array_data(i).model_switch_to_0(j);
                slope = array_data(i).model_switch_to_0_strength(j);
                plot([T-p.strength_window, T+p.strength_window], [0-slope*p.strength_window, 0+slope*p.strength_window], 'r-')
            end
            for j=1:length(array_data(i).model_switch_to_1)
                T = array_data(i).model_switch_to_1(j);
                slope = array_data(i).model_switch_to_1_strength(j);
                plot([T-p.strength_window, T+p.strength_window], [0-slope*p.strength_window, 0+slope*p.strength_window], 'k-')
            end 
        end
        if isfield(array_data(i), 'ignore_0_strength') && p.plot_strength && p.plot_ignore
            for j=1:length(array_data(i).ignore_model_switch_to_0)
                T = array_data(i).ignore_model_switch_to_0(j);
                slope = array_data(i).ignore_0_strength(j);
                plot([T-p.strength_window, T+p.strength_window], [0-slope*p.strength_window, 0+slope*p.strength_window], 'r--')
            end
            for j=1:length(array_data(i).ignore_model_switch_to_1)
                T = array_data(i).ignore_model_switch_to_1(j);
                slope = array_data(i).ignore_1_strength(j);
                plot([T-p.strength_window, T+p.strength_window], [0-slope*p.strength_window, 0+slope*p.strength_window], 'k--')
            end 
        end
    else
        if isfield(array_data(i), 'model_switch_to_0_strength') && p.plot_strength;
            for j=1:length(array_data(i).model_switch_to_0)
                T = array_data(i).model_switch_to_0(j);
                slope = array_data(i).model_switch_to_0_strength(j);
                plot([T-p.strength_window, T+p.strength_window], [0-slope, 0+slope], 'r-')
            end
            for j=1:length(array_data(i).model_switch_to_1)
                T = array_data(i).model_switch_to_1(j);
                slope = array_data(i).model_switch_to_1_strength(j);
                plot([T-p.strength_window, T+p.strength_window], [0-slope, 0+slope], 'k-')
            end 
        end
        if isfield(array_data(i), 'ignore_0_strength') && p.plot_strength && p.plot_ignore
            for j=1:length(array_data(i).ignore_model_switch_to_0)
                T = array_data(i).ignore_model_switch_to_0(j);
                slope = array_data(i).ignore_0_strength(j);
                plot([T-p.strength_window, T+p.strength_window], [0-slope, 0+slope], 'r--')
            end
            for j=1:length(array_data(i).ignore_model_switch_to_1)
                T = array_data(i).ignore_model_switch_to_1(j);
                slope = array_data(i).ignore_1_strength(j);
                plot([T-p.strength_window, T+p.strength_window], [0-slope, 0+slope], 'k--')
            end 
        end
    end
    title(['Trial ' num2str(i)])
    pause()
end

