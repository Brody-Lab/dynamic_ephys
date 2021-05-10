% based on dev_script from the accumulation_model repo 
% close all; clear all;

% load example fit and dataset
dp = set_dyn_path(1);
ex_rat = 'H084';
model_fit_fn = fullfile(dp.model_fits_dir, ['fit_analytical_' ex_rat '.mat'])

fit = fit_rat_analytical(ex_rat,'data_dir',dp.data_dir,'results_dir',dp.model_fits_dir);
load(fit.datafile)
%%

% set up parameters
TN = 1;
data(TN).pokedR  = 0; % Which choice to compute posterior for
params = fit.final;
params(1) = 0;
%params(2) = 1e-4; % accumulation noise
%params(3) = 1e-4; % sensory noise
%params(4) = 1e-4; % initial noise
params(7) = 0;
p.compute_full  = 1; 
p.compute_back  = 1;
p.compute_delta = 5;
p.save_new_ref  = 0;    
p.dt            = 1e-4;
%p.check_final   = 1;
p.b             = 50;
p.da            = 1;
p.avals         = -p.b:p.da:p.b;
p.a_sign = 2 * data(TN).pokedR - 1; % added by TB
p.da_grid       = 1; % uncommented by TB
p.a_grid        = -p.b:p.da_grid:p.b + params(7); % uncommented by TB
%%p.a_grid        = [-p.b:p.da_grid:-2 -1.75:0.25:-1 0:p.da_grid:p.b] + params(7);
%%p.a_grid        = [-p.b:2:-10 -9:p.da_grid:p.b] + params(7);
%p.a_grid_s = [p.da_grid (diff(p.a_grid(1:end-1)) + diff(p.a_grid(2:end)))/2 p.da_grid];
%% da grid notes
%% 0 bin handled elsewhere
%% fix grid sizes in compute_particles (3 places)
%% set scale factors in compute_backwards_pass
p.d_dex         = 6 * p.a_sign; % abs value of delta-mode to compute
p.n             = 10000;
p.error_tolerance = 0.0001;
computet        = 0;
plott           = 0;
particlet       = 0;

if length(unique(diff(p.avals))) > 1
    error('p.avals is not a uniform grid! I cannot verify solution for a non uniform evaluation grid. Proceed with EXTREME caution!');
end
if p.da > 1
    error('p.da > 1 may lead to unstable behavior at some point. Proceed with caution!');
end

%if p.da_grid > p.da
%    error('delta mode spacing larger than evaluation spacing, will create funky results. Proceed with caution!');
%end

if params(3) < 0.5
    error('WARNING, you might be hitting the low-variance case. See documentation for details, and for advice on how to proceed.')
end
% compute foward pass distribution
tic
[cl, cr]    = make_adapted_cat_clicks(data(TN).leftbups, data(TN).rightbups, params(5), params(6));
data(TN).clicks      = [-cl  +cr];
data(TN).times       = [data(TN).leftbups data(TN).rightbups];
data(TN).times       = round(data(TN).times/p.dt)*p.dt;
data(TN).T           = round(data(TN).T/p.dt)*p.dt;
data(TN).numsteps    = round(data(TN).T/p.dt);
data(TN).dtimes      = round(data(TN).times*(1/p.dt));
[forward] = compute_full_trial(data(TN),params,p); 
computet  = computet +toc;

% Visualize forward pass
tic
forward = compute_pdf(forward,p.avals,p);
forward.plot_zero_line = true;
forward.plot_mean_line = false;
forward.title = 'Forward';
forward.right_click_marker = '|';
forward.left_click_marker = '|';
forward.left_clicks = data(TN).leftbups;
forward.right_clicks = data(TN).rightbups;
forward.right_click_y = 25;
forward.left_click_y = -25;
forward.click_height = 5;
%%

alim = Inf*p.a_sign;
fh = figure(1); 
clf; 
subplot(2,4,1);
plot_pdf(forward);

plott = plott + toc;
% compute backwards pass
if p.compute_back;
tic
p = get_grid(data(TN), forward, params,p);
[back, posterior] = compute_backwards_pass(data(TN),params,p,forward);
computet = computet + toc;

%%
bdex = find(back.a_grid == p.d_dex); % which mode to compute
back = compute_pdf(back,p.avals,p,'mode',bdex);
posterior = compute_pdf(posterior, p.avals,p,'mode',bdex);
back.s = ones(size(back.ma))./50; %% DOES NOTHING
backF = compute_pdf(back,p.avals,p,'mixture');
posteriorF = compute_pdf(posterior, p.avals,p,'mixture');
particle = compute_particles(data(TN), params, p);
%%
% visualize backwards pass - one delta solution
tic
back.plot_zero_line = true;
back.plot_mean_line = false;
back.title = sprintf('Backward \\delta, a_N = %i',p.d_dex);
back.right_click_marker = '|';
back.left_click_marker = '|';
back.left_clicks = data(TN).leftbups;
back.right_clicks = data(TN).rightbups;
back.right_click_y = 25;
back.left_click_y = -25;
back.click_height = 5;
subplot(2,4,2);
plot_pdf(back);

% visualize posterior - one delta solution
subplot(2,4,3);
posterior.plot_zero_line = true;
posterior.plot_mean_line = false;
posterior.title = sprintf('Posterior \\delta, a_N = %i',p.d_dex);
posterior.right_click_marker = '|';
posterior.left_click_marker = '|';
posterior.left_clicks = data(TN).leftbups;
posterior.right_clicks = data(TN).rightbups;
posterior.right_click_y = 25;
posterior.left_click_y = -25;
posterior.click_height = 5;
plot_pdf(posterior);
% 
% % visualize entire backwards pass
% backF.plot_zero_line = true;
% backF.plot_mean_line = false;
% backF.title = 'Backward Full';
% backF.right_click_marker = '|';
% backF.left_click_marker = '|';
% backF.left_clicks = data(TN).leftbups;
% backF.right_clicks = data(TN).rightbups;
% backF.right_click_y = 25;
% backF.left_click_y = -25;
% backF.click_height = 5;
% subplot(2,4,5);
% plot_pdf(backF);

% visualize entire posterior 
subplot(2,4,4);
plott = plott + toc;
tic
computet = computet + toc;

tic
posteriorF.plot_zero_line = true;
posteriorF.plot_mean_line = false;
if p.a_sign > 0
    posteriorF.title = ('Posterior Full, a_N = U(0,\infty)');
else
    posteriorF.title = ('Posterior Full, a_N = U(-\infty,0)');
end
posteriorF.right_click_marker = '|';
posteriorF.left_click_marker = '|';
posteriorF.left_clicks = data(TN).leftbups;
posteriorF.right_clicks = data(TN).rightbups;
posteriorF.right_click_y = 25;
posteriorF.left_click_y = -25;
posteriorF.click_height = 5;
plot_pdf(posteriorF);
plott = plott + toc;
end

% particles
tic
particlet = particlet + toc;
tic
subplot(2,4,5);
forwardP.pdf = particle.fpdf;
forwardP.avals = particle.avals;
forwardP.T = particle.T;
forwardP.plot_zero_line = true;
forwardP.plot_mean_line = false;
forwardP.title = 'Particle Forward';
forwardP.right_click_marker = '|';
forwardP.left_click_marker = '|';
forwardP.left_clicks = data(TN).leftbups;
forwardP.right_clicks = data(TN).rightbups;
forwardP.right_click_y = 25;
forwardP.left_click_y = -25;
forwardP.click_height = 5;
plot_pdf(forwardP)

subplot(2,4,6);
backP.pdf = particle.bpdf;
backP.avals = particle.avals;
backP.T = particle.T;
backP.plot_zero_line = true;
backP.plot_mean_line = false;
%backP.title = {'Particle Backward \delta'};
backP.title = sprintf('Particle Backward \\delta, a_N = %i',p.d_dex);
backP.right_click_marker = '|';
backP.left_click_marker = '|';
backP.left_clicks = data(TN).leftbups;
backP.right_clicks = data(TN).rightbups;
backP.right_click_y = 25;
backP.left_click_y = -25;
backP.click_height = 5;
plot_pdf(backP)

subplot(2,4,7);
posteriorD.pdf = particle.dpdf;
posteriorD.avals = particle.avals;
posteriorD.T = particle.T;
posteriorD.plot_zero_line = true;
posteriorD.plot_mean_line = false;
posteriorD.title = sprintf('Particle Posterior \\delta, a_N = %i',p.d_dex);
posteriorD.right_click_marker = '|';
posteriorD.left_click_marker = '|';
posteriorD.left_clicks = data(TN).leftbups;
posteriorD.right_clicks = data(TN).rightbups;
posteriorD.right_click_y = 25;
posteriorD.left_click_y = -25;
posteriorD.click_height = 5;
plot_pdf(posteriorD)


ax = subplot(2,4,8);
posteriorP.pdf = particle.ppdf;
posteriorP.avals = particle.avals;
posteriorP.T = particle.T;
posteriorP.plot_zero_line = true;
posteriorP.plot_mean_line = false;

if p.a_sign>0
    posteriorP.title = ('Posterior Particle Full, a_N = U(0,\infty)');
else
    posteriorP.title = ('Posterior Particle Full, a_N = U(-\infty,0)');
end
posteriorP.right_click_marker = '|';
posteriorP.left_click_marker = '|';
posteriorP.left_clicks = data(TN).leftbups;
posteriorP.right_clicks = data(TN).rightbups;
posteriorP.right_click_y = 25;
posteriorP.left_click_y = -25;
posteriorP.click_height = 5;
plot_pdf(posteriorP)
plott = plott + toc;
axpos = get(ax,'position');
cb = colorbar
drawnow
set(ax,'position',axpos)
cb.Position = cb.Position + [0 0 0 -.15]
title(cb,'P(a)')
%%
allax = findall(fh,'type','axes');
set(allax,'Box','off');
ylim(allax,[-40 40]);
end_color = dp.model_color;
colormap(colormapLinear(end_color).^2)
for aa = 1:length(allax)
    caxis(allax(aa),[0 .1])
    view(allax(aa),[90 90])
end

fht = 3 * 2.5;
fw  = 4/3*2.5 * 3.75;

set(fh, 'position', [0 10 fw fht],'papersize',[fw fht],'paperpositionmode','auto')
if p.a_sign > 0
    print(fullfile(dp.fig_dir, 'dist_check_supp_right'), '-dsvg','-painters')
else
    print(fullfile(dp.fig_dir, 'dist_check_supp_left'), '-dsvg','-painters')
end

