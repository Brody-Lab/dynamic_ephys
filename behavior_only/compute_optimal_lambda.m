function [fits] = compute_optimal_lambda(fits,data,p)

% modify data to make pokedR reflect choice
for i=1:length(data);
    correct_answer = (data(i).pokedR & data(i).hit) | (~data(i).pokedR & ~data(i).hit);
    data(i).ratPokedR = data(i).pokedR;
    data(i).pokedR = correct_answer;
end

params = fits.final(2:end);
lambda_initial = fits.final(1);

% add 0 lapse rate for fits without lapse
%if length(params) == 7
%    params = [params 0];
%end

disp('Ready for optimization step')

lower_bound = -40;
upper_bound = 40;
opts.compute_full = false;
opts.dt = .001;
        
[lambda_opt, f_opt, exitflag_opt, output_opt, ~, grad_opt, hessian_opt] = ...
    fmincon(@(lambda) compute_LL(data, [lambda params], opts), ...
    lambda_initial, [], [], [], [], lower_bound,  upper_bound, [],...
    optimset('Display', 'iter-detailed'));

fits.opt_lambda = lambda_opt;
fits.opt_hessian = hessian_opt;

disp('Optimization of lambda given other parameters complete')




