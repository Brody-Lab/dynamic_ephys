function [model] = compute_model_mean(rat, sessid, varargin)

    p = inputParser;
    addParameter(p, 'fake_params', [])
    addParameter(p, 'fake_choices', [])
    addParameter(p, 'which_trials', [])
    addParameter(p, 'dt', 1e-3)
    addParameter(p, 'avals', [-10:.1:10])
    parse(p, varargin{:})
    p = p.Results;
    
    disp('computing model responses')
    
    dp = set_dyn_path;
    do_save = 1;
    save_dir = dp.model_mean_dir;
    save_mn_path = fullfile(save_dir,['model_mean_' num2str(sessid) '.mat']);
    save_posterior_path = fullfile(save_dir, ['model_posterior_' num2str(sessid) '.mat']);

    % load data
    S = load_session_data(sessid);
    data = format_data(S);
    
    if ~isempty(p.which_trials)
        data = data(p.which_trials);
    end
    
    if ~isempty(p.fake_choices)
        for tt = 1:length(data)
            data(tt).pokedR = p.fake_choices(tt);
        end
    end
    % load params
    
    if ~isempty(p.fake_params)
        params  = p.fake_params;
        do_save = 0;
    else
        fit_file = fullfile(dp.model_fits_dir, ['fit_analytical_' rat '.mat']);
        if ~exist(fit_file, 'file')
            fit = fit_rat_analytical(rat, 'data_dir', dp.behav_data_dir, ...
                'results_dir', dp.model_fits_dir);
        else
            f = load(fit_file,'fit');
            fit = f.fit;
        end
        params = fit.final;
    end
    % compute model output for each trial
   
    [model,model_p] = accumulation_model(data,params,'posterior_only',1,...
         'error_tolerance', 1e-4, 'compute_dist', 1, 'da_grid', 0.1,...
        'avals', p.avals, 'just_pdf', 1, 'eval_dt', p.dt);
   
    % save
    model_p.sessid = sessid;
    model_p.ratname = rat;
    s = whos('model');


    % save separate file with just mean because its more compact and sufficient for some analyses
    for k=1:length(model)
        model_mean(k).mean  = model(k).posterior.mean;
        model_mean(k).T     = model(k).posterior.T;
    end
    if do_save
        if s.bytes > 2e9
            save(save_posterior_path, '-v7.3', 'model','model_p')
        else
            save(save_posterior_path, 'model','model_p')
        end
        save(save_mn_path, 'model_mean','model_p')
    end    
end

