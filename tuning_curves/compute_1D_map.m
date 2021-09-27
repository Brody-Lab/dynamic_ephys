function [f1] = compute_1D_map(array_data, vec_data,p)

% need to include a time lag?

da = 0.5;
edges = -20:da:20;
centers = -19.95:da:19.95;
model_counts = zeros(1,length(centers));
spike_counts = zeros(1,length(centers));
for i=1:length(array_data)
    [n] = histcounts(array_data(i).model_mean,edges);
    model_counts = model_counts + n;

    spike_vec = array_data(i).spikes(array_data(i).spikes>0 & array_data(i).spikes<=(array_data(i).stim_end));
    spike_vec = round(spike_vec.*(1/1e-3));
    spike_vec(spike_vec == 0) = [];
    spike_vec(spike_vec > length(array_data(i).model_mean)) = [];
    spike_vec = array_data(i).model_mean(spike_vec);
    [n] = histcounts(spike_vec,edges);
    spike_counts = spike_counts + n;
end

spike_counts(spike_counts == 1) = 0;

figure; hold on;
%%model_counts = model_counts./50;
map = spike_counts./model_counts;
CI = 1.96.*sqrt((map.*(1-map))./(model_counts));
%normval = nansum(map);
%map = map./normval;
%CI = CI./normval;
shadedErrorBar(centers, map,CI,'b')
plot(centers, map, 'k')
x = mean(map(~isnan(map) & map ~=0));
plot([-20 20], [x x],'r')
ylabel('Spike Prob', 'fontsize',16);
xlabel('Model Value (a)', 'fontsize',16)
set(gca, 'fontsize',16)
ylim([0 inf])
temp = centers(find(~isnan(map)));
xlim([min(temp), max(temp)])

f1 = gcf;
