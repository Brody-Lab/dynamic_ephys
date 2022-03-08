function modelGoR = compute_model_behavior(ratname, dp, rng_in)
if nargin < 3
    rng_in = 1;
end

fit = fit_rat_analytical(ratname, 'data_dir',dp.data_dir, 'results_dir',dp.model_fits_dir);
rng(rng_in)
modelGoR = fit.pr > rand(size(fit.pr));
