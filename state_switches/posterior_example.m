sessid  = 508439; %original figure used this session
cellid  = 16905;
rat = 'H037';
% %% use to check H037 converges now
% for ff = 1:10
%     fit{ff} = fit_rat_analytical(rat,'data_dir', dp.behav_data_dir, 'reload',0);
% end

%%
dp      = set_dyn_path;

cellid  = 18181;
[sessid, rat]   = bdata('select sessid, ratname from cells where cellid={S}',cellid);
rat             = rat{1};
S               = load_session_data(sessid);
data_full       = format_data(S);
violix          = find(S.pd{1}.violations);


%%
align_strs = dyn_align_LUT;
cin_align_ind   = find(ismember(align_strs,'stimstart-cout-mask'));

%%
fit = fit_rat_analytical(rat,'results_dir',dp.model_fits_dir);
params = fit.final;
%%
d = dyn_cell_packager(cellid);
%%

% compute model output for each trial
p_in.error_tolerance    = 1e-4;
p_in.compute_dist       = 1;
p_in.da_grid            = 0.1;
p_in.da                 = 0.1;
p_in.avals              = -10:p_in.da:10;
p_in.just_pdf           = 1;
p_in.return_backwards   = 1;

nswitches = cellfun(@length,{data_full.genSwitchTimes})';
cellpref = 'l';
pd = S.pd{1};
good = pd.hits & pd.sides==cellpref & ~pd.violations;
which_trials = find(nswitches > 1 & pd.samples > 1 & good);
ii = 1
%egtrial = 51;
for tt = 2
    which_trial = which_trials(tt)
data = data_full(which_trial);
[model1,p]          = accumulation_model(data, params, 'forward', p_in);
model1(ii).backwards = compute_pdf(model1(ii).backwards,p.avals,p,'mixture');
% plot left and right clicks
model1(ii).forward.title = 'Forward';
model1(ii).backwards.title = 'Backward';

switch data(ii).pokedR
    case 0
        chosen = 'left';
        unchosen = 'right';
    case 1
        chosen = 'right';
        unchosen = 'left';
end


%%
SD  = get_sessdata(sessid);
peh = SD.peh{1};

D = model1(ii).posterior;

D.title = 'Posterior';
D.left_clicks = data(ii).leftbups;
D.right_clicks = data(ii).rightbups;
D.click_lim = 5.5;
D.click_size = 12;
D.left_click_y = 4.75;
D.right_click_y = 5.75;
D.left_click_marker = '|';
D.right_click_marker = '|';
ms = D.T(find(diff(D.mean>0)));
D.model_switches = ms(1:end);
D.state_switches = data.genSwitchTimes;
D.model_switch_y = -3.75;
D.state_switch_y = -4.75;

fig = figure(1); clf
D.plot_mean_line = 0
plot_pdf(D)

pbaspect([2.5 1 1])
ylim([-5.75 5.75 ])
fig.Position = [5 5 6 3]
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 3 3];
fig.PaperPositionMode = 'Auto';
fig.PaperSize = [3 3];
end_color = dp.model_color;

colormap(colormapLinear(end_color).^2)

%colormap(flipud(gray))
box off
ax = gca
%ax.YColor = 'w'
ax.TickDir = 'out'
caxis([0 1])
cb = colorbar
pos = ax.Position;
xlim([-0.05 max(xlim)])
cb.Position = cb.Position+[.075 0 -.01 -.275];
title(cb,{ 'posterior' 'p(a|\theta,choice)'})
ax.Position = pos;
xlabel('accumulated evidence (a)')
title('')
xlabel('time from stim onset (s)')

fig_name = fullfile(dp.fig_dir, ['posterior_' chosen '_choice']);
print(fig, fig_name,'-dsvg','-painters')
%
ts = bdata('select ts from spktimes where cellid={S}',cellid);
ts = ts{1};

cin_ts = peh(which_trial).states.cpoke1(1);
cend_ts = cin_ts + data.T;
this_ts = ts(ts > (cin_ts -5) & ts < (cend_ts + 5))';
this_ts = this_ts - cin_ts;

%

%%
% plot spikes from same trial
fh2 = figure(2); clf
fh2.Position = [5 9 6 3];
ax2 = axes;
ax2.Position = [pos(1) .25 pos(3) .2];
hold(ax2,'on');
state_switches = data.genSwitchTimes;

    
xlabel(ax2,'time from stim onset (s)');
ax2.TickDir = 'out';
box(ax2,'off');
ax2.YColor = 'w';
xlim(get(ax,'xlim'));


left_color  = [48 127 255]./255;
right_color = [0 140 54]./255 ;
plot([D.left_clicks; D.left_clicks]', 1-[0;1],  '-', ...
    'color',left_color)

plot([D.right_clicks; D.right_clicks]', 2-[0;1],  '-', ...
    'color',right_color)

plot([D.state_switches; D.state_switches],3-[0;1],...
        '-','color', 'k','linewidth',2)
plot([D.model_switches; D.model_switches],4-[0;1],...
        '-','color', 'r','linewidth',2)
    
plot(ax2,[this_ts; this_ts],5-[0;1],'k');


ylim(ax2,[0 5])


ax3 = axes;
ax3.Position = [pos(1) .55 pos(3) .25];
t = d.frate_t{cin_align_ind};
fr = d.frate{cin_align_ind}(which_trial,:);


xlim(ax3,get(ax2,'xlim'))
linkaxes([ax,ax2,ax3],'x')
ax3.XTickLabels = [];
ax3.TickDir = 'out';
krn_width = 0.1; 
bin_size = 0.005; 
dx=ceil(5*krn_width/bin_size);
krn=normpdf(-dx:dx,0,krn_width/bin_size);
krn(1:dx)=0;
krn=(krn)/sum(krn)/bin_size;
[y x] = spike_filter(cin_ts, ts, krn, 'pre', 1, 'post', 3,...
    'kernel_bin_size', bin_size, 'normalize_krn',1);
%plot(ax3,x,y)
%plot(ax3,t, fr,'k')
hold(ax3,'on')


end_state = data.genEndState;
odd_nstates = mod(length(state_switches)+1,2);

blocks = [0 D.model_switches cend_ts];
for ss = 1:length(blocks)-1
    this_t = x >= blocks(ss) & x <= blocks(ss+1);
    a = find(this_t,1,'first');
    b = find(this_t,1,'last');
    
    odd_state = end_state;
    if odd_nstates == mod(ss,2)
        this_state = end_state;
    else
        this_state = ~end_state;
    end
    if this_state
        this_color = dp.right_color;
    else
        this_color = dp.left_color;
    end
    
    patch(ax3,x([a a:b b a]),[0 y(a:b) 0 0],this_color,'edgecolor',this_color)
    plot(ax3,x(a:b),y(a:b),'color',this_color,'linewidth',1.5)
    alpha(.5)
    
    this_t_ = D.T >= blocks(ss) & D.T <= blocks(ss+1);
    a_ = find(this_t_,1,'first');
    b_ = find(this_t_,1,'last');
    plot(ax(1),D.T(a_:b_),D.mean(a_:b_),'color',this_color,'linewidth',2)

end

%

% 
% %x = t; y = fr;
plot(ax3,[D.state_switches; D.state_switches],max(ylim)-[0;10],...
        '-','color', 'k','linewidth',2)
plot(ax3,[D.model_switches; D.model_switches],max(ylim)-[0;10],...
        '-','color', 'r','linewidth',2)
box(ax3,'off')
set(ax3,'xlim',get(ax,'xlim'))

fh2_name = fullfile(dp.fig_dir,['spikes_posterior_' chosen '_choice']);
print(fh2, fh2_name,'-dsvg','-painters')

print(fig, fig_name,'-dsvg','-painters')

%%
end
%%
figure(3); clf
imagesc([d.frate{cin_align_ind}(go_r,:); d.frate{cin_align_ind}(go_l,:)])
%%
%%


model1(ii).forward.left_clicks = data(ii).leftbups;
model1(ii).forward.right_clicks = data(ii).rightbups;
model1(ii).forward.click_lim = 5.5;
model1(ii).forward.click_size = 8;
model1(ii).backwards.left_clicks = data(ii).leftbups;
model1(ii).backwards.right_clicks = data(ii).rightbups;
model1(ii).backwards.click_lim = 5.5;
model1(ii).backwards.click_size = 8;
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




