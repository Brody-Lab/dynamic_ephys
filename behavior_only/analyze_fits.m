function [F,group] = analyze_fits(group)
if isempty(strmatch(group, char('a_fits_3','a_no_sa','a_no_si','a_x_stream','synthetic','synthetic_short','dataset_2','BING','dataset_2_noL','BING_noL','dataset_2_prior','dataset_2_gaussian_prior','ephys'),'exact'))
    error(['ERROR: No set of model fits matching that group name: ' group] )
end
dp = set_dyn_path;

fn = fullfile(dp.data_dir, 'group_analysis.mat');
if exist(fn,'file')
    load(fn, 'F', 'analysis')
    return
end

disp('Here we go!')

% 0 is bing's model, 1 is my model
analytical = 1;

% Which rats to analyze
rats = {'H033','H037','H038','H039', 'H040','H043', 'H045','H058', 'H061','H065','H066','H067','H083', 'H084'};
if strcmp(group, 'BING') | strcmp(group, 'BING_noL')
    rats = {'B052', 'B053', 'B065', 'B069', 'B074', 'B083', 'B090', 'B093', 'B097', 'B102', 'B103', 'B104', 'B105', 'B106', 'B111', 'B112', 'B113', 'B115'};
end
if strcmp(group, 'ephys')
    rats = dp.ratlist;%rats = {'H037','H066','H084','H129', 'H191'};%, 'H176'};
end

%if strcmp(group, 'synthetic');
%rats = {'1','2','3','4','5','6','7','8','9'};
%end

param_names = { 'lambda',     'sa',     'ss',     'si',     'B',      'phi',    'tau',    'bias'};%,   'lapse'};


% Iterate over rats and either load or compute analysis
for i=1:length(rats)
    
    if strcmp(group, 'synthetic') | strcmp(group, 'synthetic_short')
        ratname = ['dataset_' rats{i}];
    else
        ratname = rats{i};
    end
    disp(ratname)
    
    load(fullfile(dp.behav_data_dir, ratname),'data')
    % Load Analysis for this rat
    fit = [];
    [fit, ~, fname] = fit_rat_analytical(ratname, 'results_dir', dp.model_fits_dir, ...
        'data_dir', dp.behav_data_dir);
    

    
    if ~isfield(fit,'fpt')
        % Analysis doesn't exist, lets compute it
        %try

            
            % Compute error bars on model fits
            fit = compute_parameter_uncertainty(fit,dp.model_fits_dir);
            
            if ~isreal(fit.se)
                disp('WARNING, some of the standard errors were complex!')
                fit.se = real(fit.se);
            end
            %            if 0
            if ~strcmp(group,'synthetic') & ~strcmp(group, 'synthetic_short')
                % Compute click mislocalization probability
                disp('Computing effective noise')
                fit = compute_effective_noise(fit,data);
                
                disp('Computing optimal lambda')
                fit = compute_optimal_lambda(fit,data,dp)
                
            end
            %            end
            
            % likelihood per trial
            fit.fpt = exp(-fit.f/length(data));
            fit.nt = length(data);
            
            % load excess click rate analysis
            try
                if ~(strcmp(group, 'BING') |strcmp(group, 'BING_noL') )
                    disp('Analysis complete, trying to load excess click rate')
                    load(['../check_rats/rat_data/' ratname '/excess1']);
                    fit.clicks = clicks;
                    fit.timescaleRatFit = clicks.bootRat.b;
                end
            catch
                % compute it?
                disp('Couldnt load excess click rate analysis, try checking ratname');
            end
            
            % Save analysis
            if analytical
                if exist('par','var') & ~isfield(fit, 'par') 
                    fit.par = par;
                end
                save(fname, 'fit');
            else
                save(fname, 'fit');
            end
            disp('Analysis complete, excess click rate complete, everything saved')
            
%         catch me
%             disp(me)
%             disp('Couldnt load analysis or model fit, try checking ratname')
%         end
%         
    end
    
    
    
    F{i} = fit;
    disp('Found analysis')
    
end
try
    analysis = compute_group_analysis(F);
catch
end
%%%%%%%%%%%%%%%%%%%%%%% SAVE DATA SUMMARY

save(fn, 'F','analysis');

%%%%%%%%%%%%%%%%%%%%%%% PLOTTING SECTION
%keyboard

% Analysis
%plot_optimal_lambda(F,group);
%plot_timescale_comparison(F,group);
%plot_noise_lambda_tradeoff(F,group);
%plot_parameter_comparisons(F,group)

% Older analysis
%%plot_noise_comparisons(F,group, compare_group);
%%plot_parameter_recovery(F,compare_group, group);



