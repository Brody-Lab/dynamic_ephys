function [fits] = compute_parameter_uncertainty(fits, hess_dir)

%if length(fit.final) == 8 | length(fit.final) == 7 | length(fit.final) == 6
if length(fits.final) ~= 9
    try
        autodiff_hessian = [];
        load(fullfile(hess_dir, ['/julia_hessian_' fits.rat]),'autodiff_hessian');
        fits.auto_hessian = autodiff_hessian;
        fits.auto_se = sqrt(diag(inv(autodiff_hessian)));
        %fit.matlab_se = sqrt(diag(inv(hessian)));
        fits.se = fits.auto_se;

        if cond(fits.auto_hessian) > 1e12
            disp('Matrix singular to working precision, so using matlab standard errors')
            fits.se = fits.matlab_se; 
        end
        disp('Incorporated julia hessian')
    catch
        error('Couldnt find julia hessian file')
    end
else
    h = hessian;
    h(5,:) = [];
    h(:,5) = [];
    
    se = sqrt(diag(inv(h)));
    fits.se(1:4) = se(1:4);
    fits.se(5)   = NaN;
    fits.se(6:9) = se(5:end);
    disp('Need to figure out stick bounds!!!!!!')
end
