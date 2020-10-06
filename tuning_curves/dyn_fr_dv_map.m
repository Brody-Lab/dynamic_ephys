 function results = dyn_fr_dv_map(cellids, t0s, n_dv_bins, ops, varargin)
% function fr_dv_map(cellids, t0s, n_dv_bins, varargin)
%
% Function that performs analysis of firing rate's relationship to
% behavioral decision variable. It generates two figures: 1) Firing rate versus time
% parameterized on DV value in color. 2) The FR-DV relationship averaged
% across specified times
%
% Input:
%       cellids:    vector of cellids to analyze
%       t0s:        vector of times to examine. These will be times into the
%                       model's evolution relative to alignment.
%       n_dv_bins:  number of bins to use for decision variable. Bounds are
%                   included as outside bins, so this must be >=3.
%   
%
% Output:
%       NONE
%
% 
% 
% %%% Default Parameters (change with varargin)

p = inputParser;
addParameter(p, 'lag', .2) % click lag: 0.2
addParameter(p, 'time_edges', []) % first and last time over which to calculate average; if empty use all t0s
addParameter(p, 'var_weight', true); % weight mean by inverse variance if true
addParameter(p, 'frbins', 0:0.1:50);  % bins for FR in joint distribution
addParameter(p, 'dt',  0.01); % model dt
addParameter(p, 'alignment',  'stimstart-cout-mask'); % string matching one of align_strs in dyn_align_LUT
addParameter(p, 'direction',  'backward'); % 'forward' or 'backward'
addParameter(p, 'krn_width',  []); % forces neural data to have this kernel width; if empty, uses whatever is in data
addParameter(p, 'fr_dt',  []); % forces neural data to have this bin size; if empty, uses whatever is in data
addParameter(p, 'krn_type',  'halfgauss'); % forces neural data to have this kernel type
addParameter(p, 'norm_type',  'div'); % type of firing rate normalization; 'div' is divisive, 'z' is z-score
addParameter(p, 'flip_bool',  '~data.prefsideright{align_ind}'); % boolean choice for each cellid of whether to flip sign of DV (string evaluated in workspace); default makes pref direction the higher DV value
addParameter(p, 'force_frdv',  0); % force extraction from database rather than using saved file
addParameter(p, 'force_bin',  0); 
addParameter(p, 'force_dv',  0); 
addParameter(p, 'save_path',  '');  
addParameter(p, 'save_filename',  'unknown'); % a name to specify this particular map (perhaps the group of cells used)
addParameter(p, 'save_map',  true); % whether to save the map as one big file
addParameter(p, 'axes_in',  []); % axis to use for time plot; if empty, create new figure
addParameter(p, 'axes_in_ta',  []); % axis to use for time-average plot; if empty, create new figure
addParameter(p, 'ta_color',  'k'); % time-average axis color
addParameter(p, 'ta_offset',  0); % offset in x-axis data points for time-average plot
addParameter(p, 'dv_bin_edges',  []); % if empty, makes n_dv_bins quantilized bins; otherwise, override with this 
addParameter(p, 'norm_ta',  true); % whether to "normalize" time average between zero and one
addParameter(p, 'ta_marker',  ''); % marker type for time average plot
addParameter(p, 'connect_bound',  true); % whether to connect the bound bin on plot to non-bound DVs with dotted line
addParameter(p, 'plot_errorbars',  true); % plot error bars on time plot
addParameter(p, 'ta_use_colors',  true); % whether to use different colors for different ta points     
addParameter(p, 'plot_limit',  100); % most extreme DV value to plot (usually make this correspond to lowest bound for any rat)
addParameter(p, 'Markersize',  2); % markersize of plot
addParameter(p, 'bound',  inf); % bound for slope-fitting purposes
addParameter(p, 'n_iter',  1); % number of iterations for refining estimate of DV
addParameter(p, 'param_scale_num',  1); % parameter number to scale
addParameter(p, 'param_scale_factor',  1); % multiplicative factor of that parameter
parse(p, varargin{:});
%bootstrap = false;                      % if true, sample with replacement across trials; used for bootstrap stats
struct2vars(p.Results);

if isempty(save_path)  %#ok<NODEF>
    dp          = set_dyn_path;
    save_path   = dp.celldat_dir;
    datadir     = fullfile(save_path);
end

if n_dv_bins<3
    error('You must have at least three DV bins. 1 for each bound and 1 for non-bound values.')
end

% convert time edges to bin numbers
if isempty(time_edges)
    time_bins = 1:numel(t0s);
elseif length(time_edges)~=2
    error('time_edges must be empty of a vector of length 2.')
else
    [~,first_bin]   = min(abs(t0s-time_edges(1)));
    [~,last_bin]    = min(abs(t0s-time_edges(2)));
    time_bins       = first_bin:last_bin;
end

if ~isempty(dv_bin_edges)
    n_dv_bins = length(dv_bin_edges) + 1;
end


if ~force_frdv && exist([save_path filesep save_filename '_frdvmapdata.mat'], 'file')
    load([save_path filesep save_filename '_frdvmapdata.mat']);
else
    % setup space
    n_cells             = numel(cellids);
    fga_cell            = nan(numel(t0s),n_dv_bins,n_cells);
    fvga_cell           = nan(numel(t0s),n_dv_bins,n_cells);
    fga_cell_residual   = nan(numel(t0s),n_dv_bins,n_cells);
    fga_cell_nresidual  = nan(numel(t0s),n_dv_bins,n_cells);
    fga_ta_cell         = nan(n_dv_bins,n_cells);
    fga_nta_cell        = nan(n_dv_bins,n_cells);
    fga_ta_unorm_cell   = nan(n_dv_bins,n_cells);
    fga_std_ta_cell     = nan(n_dv_bins,n_cells);
    fga_std_nta_cell    = nan(n_dv_bins,n_cells);
    fga_aa_cell         = nan(numel(t0s),n_cells);
    x_cell              = nan(n_dv_bins,n_cells);
    fr_modulation       = nan(n_cells,1);
    frm_time            = nan(numel(t0s), n_cells);
    flipdex             = zeros(n_cells,1);
    tuning_cell         = nan(numel(t0s), n_dv_bins,n_cells);
    tuning_fr           = nan(numel(t0s), n_cells);
    tuning_fa           = nan(n_dv_bins,n_cells);
    rank1_variance      = nan(n_cells,1);
   
    % loop through cells
    for ci=1:n_cells
        fprintf('Joint FR-DV analysis on cell %d of %d \n', ci, numel(cellids));
        try
%        rat = bdata('select ratname from cells where cellid="{S}"',cellids(ci));

        [x, frbins, Pjoints, fr_given_as, fr_var_given_as] = ...
            dyn_compile_binned_database(cellids(ci), t0s, n_dv_bins, ops,...
            'lag', lag, 'krn_width', krn_width, 'fr_dt', fr_dt, 'dt', dt, 'alignment', alignment,...
            'direction', direction, 'frbins', frbins, 'krn_type', krn_type, 'norm_type', norm_type,...
            'n_iter',n_iter, 'param_scale_num', param_scale_num, ...
            'param_scale_factor', param_scale_factor, 'force_bin', force_bin,'force_dv',force_dv,...
            'datadir',datadir);%, 'bootstrap', bootstrap);


        % Do we need to flip cell?
        fga_ta_temp = nanmean(fr_given_as(time_bins,:)- nanmean(fr_given_as(time_bins,:),2),1);           
        flip_cond  = mean(fga_ta_temp(x>0)) < mean(fga_ta_temp(x<0));
        if flip_cond
            fr_given_as     = flipdim(fr_given_as,2);
            fr_var_given_as = flipdim(fr_var_given_as,2);
            flipdex(ci)     = 1;
        end
        
        % Compute fga, and fga_residual 
        fga_cell(:,:,ci)    = fr_given_as;
        fga_cell_residual(:,:,ci) = fr_given_as - nanmean(fr_given_as,2);
        fga_aa_cell(:,ci)   = nanmean(fr_given_as,2);
        fvga_cell(:,:,ci)   = fr_var_given_as;
        x_cell(:,ci)        = x;
        
        % normalize fga_residual for each time bin
        for t = 1:size(fga_cell_residual,1);
            fga_cell_nresidual(t,:,ci) = fga_cell_residual(t,:,ci);
            fga_cell_nresidual(t,:,ci) = fga_cell_nresidual(t,:,ci) - min(fga_cell_nresidual(t,:,ci));
            fga_cell_nresidual(t,:,ci) = fga_cell_nresidual(t,:,ci)./max(fga_cell_nresidual(t,:,ci));
        end
        frm_time(:,ci) = max(fga_cell_residual(:,:,ci)') - min(fga_cell_residual(:,:,ci)'); 
        fga_nta     = nanmean(fga_cell_nresidual(time_bins(frm_time(:,ci) > ops.min_fr_mod),:,ci),1);
        fga_std_nta = nanstderr(fga_cell_nresidual(time_bins(frm_time(:,ci) > ops.min_fr_mod),:,ci),1);
        min_val = min(fga_nta);        
        fga_nta = fga_nta - min_val;
        max_val = max(fga_nta);
        fga_nta = fga_nta./max_val;
        fga_std_nta = fga_std_nta ./max_val;

        % Find 1D tuning curve for each cell
        fga_ta      = nanmean(fr_given_as(time_bins,:)- nanmean(fr_given_as(time_bins,:),2),1);
        fga_std_ta  = nanstderr(fr_given_as(time_bins,:)-nanmean(fr_given_as(time_bins,:),2),1);
        fga_ta_unorm = fga_ta;
        min_val     = min(fga_ta);
        fga_ta      = fga_ta - min_val;
        max_val     = max(fga_ta);
        fga_ta      = fga_ta ./ max_val;
        fga_std_ta  = fga_std_ta ./ max_val;
        fga_ta_cell(:,ci)       = fga_ta;
        fga_nta_cell(:,ci)      = fga_nta;
        fga_std_ta_cell(:,ci)   = fga_std_ta;
        fga_std_nta_cell(:,ci)   = fga_std_nta;
        fga_ta_unorm_cell(:,ci) = fga_ta_unorm;
        fr_modulation(ci) = max(fga_ta_unorm) - min(fga_ta_unorm);

        % fit sigmoid to 1D tuning curve - using direct averaging
        try
            warning('off','all')
            [betas,resid,jacob,sigma,mse] = nlinfit(x,fga_ta,@dyn_sig,[0, 1/3]);
            warning('on','all')
            betas_cell(ci,:)    = betas;
            delta               = sqrt(diag(sigma)) * tinv(1-0.05/2,sum(~isnan(x))-4);
            sigmas_cell(ci,:)   = delta;
            slope_cell(ci)      = betas(2)/4;
        catch
            betas_cell(ci,:)    = nan;
            sigmas_cell(ci,:)   = nan;
            slope_cell(ci)      = nan;
        end
        % fit sigmoid to 1D tuning curve - using normalized / time bin
        try
            warning('off','all')
            [betas,resid,jacob,sigma,mse] = nlinfit(x,fga_nta,@dyn_sig,[0, 1/3]);
            warning('on','all')
            nbetas_cell(ci,:)    = betas;
            ndelta               = sqrt(diag(sigma)) * tinv(1-0.05/2,sum(~isnan(x))-4);
            nsigmas_cell(ci,:)   = delta;
            nslope_cell(ci)      = betas(2)/4;
        catch
            nbetas_cell(ci,:)    = nan;
            nsigmas_cell(ci,:)   = nan;
            nslope_cell(ci)      = nan;
        end
        % fit sigmoid to 1D tuning curve for each normalized time bin
        if ops.fit_time_sigmoid
        warning('off','all')
        for t = 1:size(fga_cell_nresidual,1);
            try
                [betas,resid,jacob,sigma,mse] = nlinfit(x,fga_cell_nresidual(t,:,ci),@dyn_sig,[0, 1/3]);
                nbetas_time(ci,t,:)    = betas;
                ndelta                      = sqrt(diag(sigma)) * tinv(1-0.05/2,sum(~isnan(x))-4);
                nsigmas_time(ci,t,:)   = delta;
                nslope_time(ci,t)      = betas(2)/4;
            catch
                nbetas_time(ci,t,:)    = nan;
                nsigmas_time(ci,t,:)   = nan;
                nslope_time(ci,t)      = nan;
            end
        end
        warning('on','all')
        end

        % SVD analysis
        [u,s,v] = svd(fga_cell_residual(:,:,ci));
        s2 = s;
        s2(2:end) = 0;
        tuning_cell(:,:,ci) = u*s2*v';
        alpha   = 1/range(v(:,1));
        beta    = s2(1,1)/alpha;
        tuning_fr(:,ci) = u(:,1)*beta;
        tuning_fa(:,ci) = v(:,1)*alpha;
        s_squared = diag(s).^2;
        rank1_variance(ci) = s_squared(1)./sum(s_squared);
        tuning_alpha(ci) = alpha;
        tuning_beta(ci) = beta;  
        if mean(tuning_fr(:,ci)) < 0
        tuning_fr(:,ci) = -u(:,1)*beta;
        tuning_fa(:,ci) = -v(:,1)*alpha;
        end 
        % fit sigmoid to rank 1 tuning curve
        try
            warning('off','all')
            svdta = tuning_fa(:,ci)';
            svdta = svdta - min(svdta);
            svdta = svdta./max(svdta); 
            [betas,resid,jacob,sigma,mse] = nlinfit(x,svdta,@dyn_sig,[0, 1/3]);
            warning('on','all')
            svd_betas_cell(ci,:)    = betas;
            svd_delta               = sqrt(diag(sigma)) * tinv(1-0.05/2,sum(~isnan(x))-4);
            svd_sigmas_cell(ci,:)   = svd_delta;
            svd_slope_cell(ci)      = betas(2)/4;
        catch
            svd_betas_cell(ci,:)    = nan;
            svd_sigmas_cell(ci,:)   = nan;
            svd_slope_cell(ci)      = nan;
        end

        catch me
            disp('something crashed in dyn_fr_dv_map')
            keyboard
        end    
    end
    
    results.fvga_cell           = fvga_cell;
    results.fga_cell            = fga_cell;
    results.fga_cell_residual   = fga_cell_residual;
    results.fga_cell_nresidual  = fga_cell_nresidual;
    results.fga_ta_cell         = fga_ta_cell;
    results.fga_nta_cell        = fga_nta_cell;
    results.fga_std_ta_cell     = fga_std_ta_cell;
    results.fga_std_nta_cell    = fga_std_nta_cell;
    results.fga_ta_unorm_cell   = fga_ta_unorm_cell;
    results.fga_aa_cell         = fga_aa_cell;
    results.betas_cell          = betas_cell;
    results.sigmas_cell         = sigmas_cell; 
    results.slope_cell          = slope_cell;
    results.nbetas_cell         = nbetas_cell;
    results.nsigmas_cell        = nsigmas_cell; 
    results.nslope_cell         = nslope_cell;
    if ops.fit_time_sigmoid
    results.nbetas_time         = nbetas_time;
    results.nsigmas_time        = nsigmas_time; 
    results.nslope_time         = nslope_time;
    end
    results.fr_modulation       = fr_modulation;
    results.frm_time            = frm_time;
    results.flipdex             = flipdex;
    results.x_cell              = x_cell;
    results.tuning_cell         = tuning_cell;
    results.tuning_fr           = tuning_fr;
    results.tuning_fa           = tuning_fa;
    results.rank1_variance      = rank1_variance;
    results.svd_betas_cell      = svd_betas_cell;
    results.svd_sigmas_cell     = svd_sigmas_cell; 
    results.svd_slope_cell      = svd_slope_cell;
    results.tuning_alpha        = tuning_alpha;
    results.tuning_beta         = tuning_beta;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finished iterating across cells, now population analysis
ops.var_weight        = var_weight;
ops.cell_dex          = results.fr_modulation > ops.min_fr_mod;
results.cell_dex    = ops.cell_dex;
results.cell_num    = find(ops.cell_dex);
results.t0s         = t0s;
results.dv_axis     = mean(results.x_cell,2);
results.time_bins   = time_bins;
results.cellid      = cellids;

try
    results             = population_analysis(results,ops);
catch me
    disp(me.message);
end
% plot population summary figures, only plots good cells
try
    plot_population(results);
catch me
    disp(me.message);
end

numgood = sum(ops.cell_dex);
fracgood = 100*(numgood/numel(cellids));
disp(['Had ' num2str(numgood) ' good cells out of ' num2str(numel(cellids)) ' : ' num2str(fracgood) '%'])

% Save
if save_map,
    if ~exist(save_path, 'dir'), mkdir(save_path); end;
    save([save_path filesep save_filename '_frdvmapdata.mat'], 'results');
end



