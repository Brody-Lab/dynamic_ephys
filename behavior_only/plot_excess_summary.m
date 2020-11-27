clear
dp = set_dyn_path;
gapath = fullfile(dp.data_dir, 'group_analysis.mat');
load(gapath)

for ff=1:length(F)
    tic
    fprintf('working on rat %i (of %i)', ff, length(F))
    close all
    if isfield(F{ff}, 'clicks')
        plot_excess_rates(F{ff}.clicks)
        print(fullfile(dp.fig_dir, [F{ff}.rat '_excess']), '-dsvg')
    else
        load(fullfile(dp.model_fits_dir, ['fit_analytical_' F{ff}.rat '.mat']));
        load(fullfile(dp.data_dir, F{ff}.rat))

        p = struct();
        clicks = compute_excess_rates_data(data,p) 
        F{ff}.clicks = clicks;
    end
    toc
end

save(gapath,'F','analysis')


