function [results] = population_analysis(results,p)

% compute only for certain cells
if isfield(p, 'cell_dex')
    fvga_cell           = results.fvga_cell(:,:,p.cell_dex);
    fga_cell            = results.fga_cell(:,:,p.cell_dex);
    if p.use_nresidual
        fga_cell_residual   = results.fga_cell_nresidual(:,:,p.cell_dex);
    else
        fga_cell_residual   = results.fga_cell_residual(:,:,p.cell_dex);
    end

    fga_ta_cell         = results.fga_ta_cell(:,p.cell_dex);
    fga_std_ta_cell     = results.fga_std_ta_cell(:,p.cell_dex);
    fga_aa_cell         = results.fga_aa_cell(:,p.cell_dex);
    betas_cell          = results.betas_cell(p.cell_dex,:);
    sigmas_cell         = results.sigmas_cell(p.cell_dex,:);
    slope_cell          = results.slope_cell(p.cell_dex);
    tuning_fa_cell      = results.tuning_fa(:,p.cell_dex);
else
    fvga_cell           = results.fvga_cell;
    fga_cell            = results.fga_cell;
    if p.use_nresidual
        fga_cell_residual   = results.fga_cell_nresidual;
    else
        fga_cell_residual   = results.fga_cell_residual;
    end
    fga_ta_cell         = results.fga_ta_cell;
    fga_std_ta_cell     = results.fga_std_ta_cell;
    fga_aa_cell         = results.fga_aa_cell;
    betas_cell          = results.betas_cell;
    sigmas_cell         = results.sigmas_cell;
    slope_cell          = results.slope_cell;
    tuning_fa_cell      = results.tuning_fa;
end
n_cells = size(fga_ta_cell,2);

% inverse variance weighted mean of cells or straight mean
if ~p.var_weight
    % mean across cells 
    ymn     = nanmean(fga_cell_residual,3); 
    %yvar    = nansum(fvga_cell./(n_cells^2),3);
    %yst     = sqrt(yvar);
    %yst2    = nanstd(fga_cell_residual,0,3) ./ sqrt(n_cells);  
        
    % average time and cells simultaneously
    data    = fga_cell_residual(results.time_bins,:,:);
    data    = shiftdim(data,1);
    data    = reshape(data,size(data,1),size(data,2)*size(data,3));
    ymn_ta  = nanmean(data,2)';
  %  yst_ta  = nanstderr(data,2)';
    yst_ta  = nanstd(data,0,2)'./ sqrt(n_cells);
   
    % cell weights
    results.w = [];
    [u,s,v] = svd(ymn);
    s2 = s;
    s2(2:end) = 0;
    tuning_population = u*s2*v';
    alpha       = 1/range(v(:,1));
    beta        = s2(1,1)/alpha;
    tuning_fr   = u(:,1)*beta;
    tuning_fa   = v(:,1)*alpha;
    s_squared   = diag(s).^2;
    rank1_variance   = s_squared(1)./sum(s_squared);
    tuning_alpha     = alpha;
    tuning_beta      = beta;  
    if mean(tuning_fr) < 0
        tuning_fr = -u(:,1)*beta;
        tuning_fa = -v(:,1)*alpha;
    end 

else
    error('you havent dealt with this case yet')
    % mean across cells: weight by variance (at each timepoint) 
    fr_bin_size = frbins(2) - frbins(1);
    var_min = 0;
    % bin variances based on resolution of FR bins
    % fvga_cell = round(fvga_cell/fr_bin_size)*fr_bin_size;
    fvga_cell = floor(fvga_cell/fr_bin_size)*fr_bin_size;
    fvga_cell(fvga_cell<=var_min) = nan;   % zero variances are data-deficient

    % fvga_cell_mn = nanmean(nanmean(fvga_cell,2),1);
    % fvga_cell = repmat(fvga_cell_mn,size(fvga_cell,1),size(fvga_cell,2));
    w_denom = nansum(1./fvga_cell,3);
    w       = (1./fvga_cell) ./ repmat(w_denom,[1 1 n_cells]);
    results.w = w;
    ymn     = nansum(fga_cell .* w,3);
    yvar    = nansum(fvga_cell .* (w.^2),3); % FIX: try using nanvar of distribution instead of fvga_cell for this
    yst     = sqrt(yvar);
    yvar2   = nansum(repmat(nanstd(fga_cell,0,3),[1 1 n_cells]) .* (w.^2),3);
    yst2    = sqrt(yvar2);  % non-propagation of variance method
    
    % mean across time, variance weighted method 
    w_ta    = w./length(results.time_bins);
    ymn_ta  = nansum(nansum(fga_cell(results.time_bins,:,:).*w_ta(results.time_bins,:,:),3),1);  % average time after across cells
    yvar_ta = nansum(nansum(fvga_cell(results.time_bins,:,:).*(w_ta(results.time_bins,:,:).^2),3),1);
    yst_ta  = sqrt(yvar_ta);
end

% Fit the sigmoid
min_val = min(ymn_ta);
ymn_ta = ymn_ta - min_val;
max_val = max(ymn_ta);
ymn_ta = ymn_ta ./ max_val;
yst_ta = yst_ta ./ max_val;
[betas,resid,jacob,sigma,mse] = nlinfit(results.dv_axis',ymn_ta,@dyn_sig,[0, 1/3]);

results.fga_population      = ymn;
results.fga_ta_population   = ymn_ta;
results.fga_std_ta_pop      = yst_ta;
results.grand_fit_betas     = betas;
results.grand_fit_sigmas    = nlparci(betas,resid,'covar',sigma);
results.grand_fit_covar     = sigma;
results.grand_fit_resid     = resid;
results.grand_slope         = betas(2)/4;
results.tuning_population   = tuning_population;
results.tuning_pop_fr           = tuning_fr;
results.tuning_pop_fa           = tuning_fa;
results.tuning_pop_r1_var       = rank1_variance;





