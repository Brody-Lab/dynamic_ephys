sessid  = 508439; %original figure used this session
cellid  = 16905;
rat = 'H037';
lag = .1;
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
ii = 1;
%egtrial = 51;
tt = 2;
which_trial = which_trials(tt)
data = data_full(which_trial);
[model1,p]          = accumulation_model(data, params, ...
    'return_backwards',1,'forward', p_in, 'compute_dist', 1);
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
sp = struct('clear_bad_strengths',1,'bad_strength',0,'fit_line',1,...
    'exclude_final',0,'final_only',0,...
    'min_pre_dur',0,'min_post_dur',0,...
    'model_smooth_wdw',100,'t_buffers',[.2 .2]);

[~, ~, array_data, vec_data] = ...
    get_switches(cellid, ...
    'which_switch','model',...
    'clear_bad_strengths', sp.clear_bad_strengths, ...
    'bad_strength', sp.bad_strength, 'fit_line', sp.fit_line,...
    'exclude_final', sp.exclude_final, 'final_only', sp.final_only,...
    'min_pre_dur',sp.min_pre_dur,'min_post_dur',sp.min_post_dur,...
    'model_smooth_wdw', sp.model_smooth_wdw,...
    't_buffers', sp.t_buffers);

tn = find([array_data.trialnum]==which_trial);
ad = array_data(tn);

%% plot accumulator distribution and clicks for this trial
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
D_ = D;
D_ = rmfield(D_,'model_switches');
D_ = rmfield(D_, 'state_switches');
D_ = rmfield(D_, 'left_clicks');
D_ = rmfield(D_, 'right_clicks');
plot_pdf(D_)

pbaspect([2.5 1 1])
ylim([-5.75 5.75 ])

fht = 2.5;
fw  = 3.75;
set(fig,'position',[5 5 fw fht],'papersize', [fw fht])
set(gca,'fontsize',dp.fsz)

% fig.Position = [5 5 6 3]
% fig.PaperUnits = 'inches';
% fig.PaperPosition = [0 0 3 3];
% fig.PaperPositionMode = 'Auto';
% fig.PaperSize = [3 3];
end_color = dp.model_color;

colormap(colormapLinear(end_color).^2)

%colormap(flipud(gray))
box off
ax = gca
%ax.YColor = 'w'
ax.TickDir = 'out'
caxis([0 1])

pos = ax.Position;
xlim([-0.05 max(xlim)])
cb = colorbar

cb.Position = cb.Position+[.0 .35 -.025 -.25];
title(cb,{ 'posterior' 'p(a|\theta,choice)'})
ax.Position = pos;
ylabel('accumulated evidence (a)')
title('')
xlabel('time from stim onset (s)')

%%
ts = bdata('select ts from spktimes where cellid={S}',cellid);
ts = ts{1};

cin_ts = peh(which_trial).states.cpoke1(1);
cend_ts = cin_ts + data.T;
this_ts = ts(ts > (cin_ts -5) & ts < (cend_ts + 5))';
this_ts = this_ts - cin_ts;

%% plot spikes and smoothed rates from same trial
fh2 = figure(2); clf
fh2.Position = [5 9 6 3];
ax2 = axes;
ax2.Position = [pos(1) .125 pos(3) .2];
hold(ax2,'on');
state_switches = data.genSwitchTimes;

D.model_switches = sort([ad.model_switch_to_0 ad.model_switch_to_1]);


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
%ax3.XTickLabels = [];
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

model_mean = movmean(D.mean,100);

end_state = data.genEndState;
odd_nstates = mod(length(state_switches)+1,2);

blocks = [0 D.model_switches cend_ts];
blocks_phys = [0 D.model_switches + lag cend_ts] ;
for ss = 1:length(blocks)-1
    this_t = x >= blocks_phys(ss) & x <= blocks_phys(ss+1);
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
    plot(ax3,x(a:b),y(a:b),'color',this_color,'linewidth',1)
    alpha(.5)
    
    this_t_ = D.T >= blocks(ss) & D.T <= blocks(ss+1);
    a_ = find(this_t_,1,'first');
    b_ = find(this_t_,1,'last');
    plot(ax(1),D.T(a_:b_),model_mean(a_:b_),...
        'color',this_color,'linewidth',1)
    
end


plot(ax3,[D.state_switches; D.state_switches],max(ylim)-[0;10],...
    '-','color', 'k','linewidth',2)
plot(ax3,[D.model_switches; D.model_switches],max(ylim)-[0;10],...
    '-','color', 'r','linewidth',2)
box(ax3,'off')
set(ax3,'xlim',get(ax,'xlim'))
ylabel(ax3,'spikes')
xlabel(ax3,'time from stim onset (s)');

ylim(ax(1),[-1 1 ].*5.75)
%%
ylim(ax3,[-5 100])
fht = 2.5;
fw  = 3.5;

axpos = get(ax3,'position')
set(fig, 'position', [5 5 fw fht])
set(fh2, 'position', [5 10 fw fht])
%%

%%
fh2_name = fullfile(dp.fig_dir,['spikes_posterior_' chosen '_choice']);
fig_name = fullfile(dp.fig_dir,['model_posterior_' chosen '_choice']);

print(fh2, fh2_name,'-dsvg','-painters')
print(fig, fig_name,'-dsvg','-painters')

