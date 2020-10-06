function [model] = compute_model_mean(rat, sessid)

    disp('computing model responses')
    
    dp = set_dyn_path;
    save_dir = dp.model_mean_dir;
    save_mn_path = fullfile(save_dir,['model_mean_' num2str(sessid) '.mat']);
    save_posterior_path = fullfile(save_dir, ['model_posterior_' num2str(sessid) '.mat']);

    % load data
    S = load_session_data(sessid);
    data = format_data(S);

    % load params
    fit_file = fullfile(dp.model_fits_dir, ['fit_analytical_' rat '.mat']);
    if ~exist(fit_file, 'file')
        fit_rat_analytical(rat, 'data_dir', dp.behav_data_dir, ...
            'results_dir', dp.model_fits_dir);
    else
       f = load(fit_file,'fit');
       fit = f.fit;
    end
    params = fit.final;

    % compute model output for each trial
    p_in.error_tolerance    = 1e-4;
    p_in.compute_dist       = 1;
    p_in.da_grid            = 0.1;
    p_in.da                 = 0.1;
    p_in.avals              = -10:p_in.da:10;
    p_in.just_pdf           = 1;
    [model,model_p] = accumulation_model(data,params,'posterior',p_in);

    % save
    model_p.sessid = sessid;
    model_p.ratname = rat;
    s = whos('model');
    if s.bytes > 2e9
        save(save_posterior_path, '-v7.3', 'model','model_p')
    else
        save(save_posterior_path, 'model','model_p')
    end

    % save separate file with just mean because its more compact and sufficient for some analyses
    for k=1:length(model)
        model_mean(k).mean  = model(k).posterior.mean;
        model_mean(k).T     = model(k).posterior.T;
    end
    save(save_mn_path, 'model_mean','model_p')
end

