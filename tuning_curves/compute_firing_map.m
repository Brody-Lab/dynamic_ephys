function [] = compute_firing_map(array_data, vec_data, p)
% need to use spike density, to control for variable trial duration

figure; hold on;
% Iterate over trials

if p.firing_map_plot_trajectories
for i=1:length(array_data)
    % plot trajectories
    tvec = 0:1e-4:array_data(i).stim_end;
    tvec = tvec(1:length(array_data(i).model_mean));
    plot(tvec, array_data(i).model_mean, 'k--')
end
end

for i=1:length(array_data) 
    % Iterate over spikes
    for j=1:length(array_data(i).spikes)
        % check if spike train is in the right window
        if array_data(i).spikes(j) >= 0 && array_data(i).spikes(j) <= array_data(i).stim_end
            timedex = array_data(i).spikes(j);
            dex = round(timedex*(1/1e-3))+1;
            if dex <= length(array_data(i).model_mean)
            adex = array_data(i).model_mean(dex);
            plot(timedex, adex, 'ko', 'markerfacecolor','k') 
            end
        end
    end
end

