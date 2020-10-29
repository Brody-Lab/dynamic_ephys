dp      = set_dyn_path;
sessid  = 508439;
S       = load_session_data(sessid);
data    = format_data(S);
rat     = S.ratname{1};
%%

fitdata = load('/Users/oroville/projects/pbups_dyn/check_rats/ephys_data/ephys_1_/model_data/H037.mat')

%%
f = fit_rat_analytical(fitdata.data);
%%
f_other.fit = f;
%%
fty = fit_rat_analytical(rat,'results_dir',dp.model_fits_dir);
%%
fit = load('/Users/oroville/Dropbox/model_fits/ANALYSIS/ephys/fit_analysis_analyticalH037.mat');
%%
params  = fit.fit.final;

% compute model output for each trial
p_in.error_tolerance    = 1e-4;
p_in.compute_dist       = 1;
p_in.da_grid            = 0.1;
p_in.da                 = 0.1;
p_in.avals              = -10:p_in.da:10;
p_in.just_pdf           = 1;
p_in.return_backwards   = 1;

ii = 1;
data(2:end) = [];
[model1,p]          = accumulation_model(data, params, 'forward', p_in);
model1(ii).backwards = compute_pdf(model1(ii).backwards,p.avals,p,'mixture');
% plot left and right clicks
model1(ii).posterior.title = 'Posterior';
model1(ii).forward.title = 'Forward';
model1(ii).backwards.title = 'Backward';
model1(ii).posterior.left_clicks = data(ii).leftbups;
model1(ii).posterior.right_clicks = data(ii).rightbups;
model1(ii).posterior.click_lim = 5.5;
model1(ii).posterior.click_size = 8;
model1(ii).forward.left_clicks = data(ii).leftbups;
model1(ii).forward.right_clicks = data(ii).rightbups;
model1(ii).forward.click_lim = 5.5;
model1(ii).forward.click_size = 8;
model1(ii).backwards.left_clicks = data(ii).leftbups;
model1(ii).backwards.right_clicks = data(ii).rightbups;
model1(ii).backwards.click_lim = 5.5;
model1(ii).backwards.click_size = 8;
%%
switch data(ii).pokedR
    case 0
        chosen = 'left';
        unchosen = 'right';
    case 1
        chosen = 'right';
        unchosen = 'left';
end
fig = figure(1); clf
plot_pdf(model1(ii).posterior)
pbaspect([3 1 1])
ylim([-6.25 6.25 ])
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 6];
fig.PaperPositionMode = 'Manual';
fig.PaperSize = [6 6];
end_color = [1 0 1]*.5;

colormap(colormapLinear(end_color).^2)
%colormap(flipud(gray))
box off
ax = gca
ax.YColor = 'w'
ax.TickDir = 'out'
caxis([0 1])
cb = colorbar
pos = ax.Position;

cb.Position = cb.Position+[.00 0 -.01 -.1];
ylabel(cb,'p(a|model)')
ax.Position = pos;
xlabel('time from stim onset (s)')
print(fig, ['posterior_' chosen '_choice'],'-dsvg')
%%
caxis([0 .25])
%%
fig = figure(3);
plot_pdf(model1(ii).backwards)
pbaspect([1 2 1])
xlim([-5.75 5.75])
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 6];
fig.PaperPositionMode = 'Manual';
fig.PaperSize = [6 6];
caxis([0 .0001])
colormap(hot)
print(['backwards_' chosen '_choice'],'-dsvg')






%%
data(ii).pokedR     = 0;
[model2,p]          =accumulation_model(data, params, 'forward', p_in);
model2(ii).backwards = compute_pdf(model2(ii).backwards,p.avals,p,'mixture');

model2(ii).posterior.title = 'Posterior';
model2(ii).forward.title = 'Forward';
model2(ii).backwards.title = 'Backward';
model2(ii).posterior.left_clicks = data(ii).leftbups;
model2(ii).posterior.right_clicks = data(ii).rightbups;
model2(ii).posterior.click_lim = 5.5;
model2(ii).posterior.click_size = 8;
model2(ii).forward.left_clicks = data(ii).leftbups;
model2(ii).forward.right_clicks = data(ii).rightbups;
model2(ii).forward.click_lim = 5.5;
model2(ii).forward.click_size = 8;
model2(ii).backwards.left_clicks = data(ii).leftbups;
model2(ii).backwards.right_clicks = data(ii).rightbups;
model2(ii).backwards.click_lim = 5.5;
model2(ii).backwards.click_size = 8;

figure(2);
plot_pdf(model2(ii).posterior)
pbaspect([1 2 1])
xlim([-5.75 5.75])
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 6];
fig.PaperPositionMode = 'Manual';
fig.PaperSize = [6 6];
print(['posterior_' unchosen '_choice'],'-dsvg')
%%
figure(4);
plot_pdf(model2(ii).backwards)
pbaspect([1 2 1])
xlim([-5.75 5.75])
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 6];
fig.PaperPositionMode = 'Manual';
fig.PaperSize = [6 6];
print(['backwards_' unchosen '_choice'],'-dsvg')


figure(5);
plot_pdf(model2(ii).forward)
pbaspect([1 2 1])
xlim([-5.75 5.75])
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 6];
fig.PaperPositionMode = 'Manual';
fig.PaperSize = [6 6];
print(['forward_' unchosen '_choice'],'-dsvg')




