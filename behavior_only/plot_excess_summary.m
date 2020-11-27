clear
dp = set_dyn_path;
gapath = fullfile(dp.data_dir, 'group_analysis.mat');
load(gapath)

for i=1:length(F)
    close all
    if isfield(F{i}, 'clicks')
        plot_excess_rates(F{i}.clicks)
        print(fullfile(dp.fig_dir, [F{i}.rat '_excess']), '-dsvg')
    else
        load(fullfile(dp.model_fits_dir, ['fit_analytical_' F{i}.rat '.mat']));
        load(fullfile(dp.data_dir, F{i}.rat))

        p = struct();
        clicks = compute_excess_rates_data(data,p) 
        F{i}.clicks = clicks;
    end
end

save(gapath,'F','analysis')


