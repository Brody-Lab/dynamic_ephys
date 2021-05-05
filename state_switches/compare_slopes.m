psg = load(['/Users/oroville/projects/pbups_dyn/code/dynamic_ephys/tuning_curves/'...
    'post_switch_gain_workspace.mat']);
%%
model_prefix = 'modelswitch_';
fn = @(fig_prefix) fullfile(dp.model_fits_dir, ['pop_' fig_prefix...
    'tuning_res.mat']);
pop_res_on = load(fn(''),'pop_res');
pop_res_on = pop_res_on.pop_res;

pop_res_ms = load(fn(model_prefix),'pop_res');
pop_res_ms = pop_res_ms.pop_res;

subtract_min = 1;
%% add population average
dx = pop_res_on(1).dv_axis;
ra_pop_on = [pop_res_on.(ra_field)]';
flip_cell = mean(ra_pop_on(:,dx>0),2) < mean(ra_pop_on(:,dx<0),2);
pop_fgta = cat(3,pop_res_on.fr_given_ta);
pop_fgta_r = cat(3,pop_res_on.fgta_resid);
ra_pop_on(flip_cell,:) = fliplr(ra_pop_on(flip_cell,:));
pop_fgta_r(:,:,flip_cell) = fliplr(pop_fgta_r(:,:,flip_cell)); 
pop_fgta(:,:,flip_cell) = fliplr(pop_fgta(:,:,flip_cell)); 
mean_pop_fgta_r = mean(pop_fgta_r,3);
mean_pop_fgta = mean(pop_fgta,3);
tempres = pop_res_on(1);
tempres.fgta_resid = mean_pop_fgta_r;
tempres.fr_given_ta = mean_pop_fgta;
tempres = compute_rank1_fgta_approx(tempres,'which_map','fr_given_ta');
ra_pop = tempres.(ra_field)';
if subtract_min
    ra_pop = ra_pop - min(ra_pop);
end
bpop_on  = fit_four_param_psycho(dx,ra_pop);

bpop_on_slope = prod(bpop_on([2 3]))/4;
ypop_on = fourParamPsychometric(bpop_on, xx);
%%
dx = pop_res_ms(1).dv_axis;
ra_pop_ms = [pop_res_ms.(ra_field)]';
flip_cell = mean(ra_pop_ms(:,dx>0),2) < mean(ra_pop_ms(:,dx<0),2);
pop_fgta = cat(3,pop_res_ms.fr_given_ta);
pop_fgta_r = cat(3,pop_res_ms.fgta_resid);
ra_pop_ms(flip_cell,:) = fliplr(ra_pop_ms(flip_cell,:));
pop_fgta_r(:,:,flip_cell) = fliplr(pop_fgta_r(:,:,flip_cell)); 
pop_fgta(:,:,flip_cell) = fliplr(pop_fgta(:,:,flip_cell)); 
mean_pop_fgta_r = mean(pop_fgta_r,3);
mean_pop_fgta = mean(pop_fgta,3);
tempms = pop_res_ms(1);
tempms.fgta_resid = mean_pop_fgta_r;
tempms.fr_given_ta = mean_pop_fgta;
tempms = compute_rank1_fgta_approx(tempms,'which_map','fr_given_ta');
ra_ms = tempms.(ra_field)';
if subtract_min
    ra_ms = ra_ms - min(ra_ms);
end
bpop_ms  = fit_four_param_psycho(dx,ra_ms);

bpop_ms_slope = prod(bpop_ms([2 3]))/4;
ypop_ms = fourParamPsychometric(bpop_ms, xx);

%%

figure(1); clf
nc = length(pop_res_on);
ncols = 8;
nrows = ceil(nc/ncols);
ra_field = 'rank1_ra_n';
fa_field = 'fga_tmn_n';

br_on = nan(nc,4);
br_ms = nan(nc,4);
bf_on = nan(nc,4);
bf_ms = nan(nc,4);
for cc = 1:nc
    %%
cellid = pop_res_on(cc).cellid;
psg_ix = find(psg.pop_cellids == cellid);

assert(cellid == pop_res_ms(cc).cellid)

diff_p = psg.diff_p(psg_ix);

dx = pop_res_on(cc).dv_axis';
xx = linspace(dx(1),dx(end),100);
ra_on = pop_res_on(cc).(ra_field);
ra_ms = pop_res_ms(cc).(ra_field);
if subtract_min
    ra_on = ra_on - min(ra_on);
    ra_ms = ra_ms - min(ra_ms);
end
br_on(cc,:)  = fit_four_param_psycho(dx,ra_on);
yy_on = fourParamPsychometric(br_on(cc,:), xx);
br_ms(cc,:)  = fit_four_param_psycho(dx,ra_ms);
yy_ms = fourParamPsychometric(br_ms(cc,:), xx);

subplot(ncols, nrows, cc);
plot(xx,yy_on,'k')
hold on
plot(xx,yy_ms,'r')
plot(dx,ra_on,'.k')
plot(dx,ra_ms,'.r')
title(cellid)
if subtract_min
    ylim([0 1])
else
    ylim([-.55 .55])
end
drawnow

if diff_p > .95
    set(gca,'xcolor','m','ycolor','m')
end
%pause(.5)
end
%%
fh = figure(2); clf
set(fh,'position',[5 5 dp.fw dp.fw]*1.5)
ax = axes;
frm = [pop_res_on.fr_mod]';

slope_on = prod(br_on(:,[2 3]),2)./4;
slope_ms = prod(br_ms(:,[2 3]),2)./4;

lims = [-1 1].*.3;
set(gca,'fontsize',dp.fsz)

plot(lims,lims,'-','color',[1 1 1].*.8)
hold on

plot([ 0 0],lims,'-','color',[1 1 1].*.8)
plot(lims,[0 0],'-','color',[1 1 1].*.8)
k = 20;
scatter(slope_on(frm<1),slope_ms(frm<1),k*frm(frm<1))
scatter(slope_on(frm>1),slope_ms(frm>1),k*frm(frm>1),'filled',...
    'markerfacecolor',[1 1 1].*.7)
scatter(slope_on(psg.diff_p>.95),slope_ms(psg.diff_p>.95),...
    k*frm(psg.diff_p>.95),'markerfacecolor','m')
%plot(slope_on(psg.diff_p>.95),slope_ms(psg.diff_p>.95),'.m')
colormap(colormapRedBlue)
xlim(lims)
ylim(lims)

xlabel('ONSET-ALIGNED slope (of rank 1 r(a))')
ylabel('SWITCH-TRIGGERED slope (of rank 1 r(a))')
box off
axis square

plot(bpop_on_slope, bpop_ms_slope, 'r+','linewidth',3,'markersize', 10,...
    'markerfacecolor',[1 1 1].*0,'markeredgecolor',[1 0 0])

%%
fh = figure(2); clf
set(fh,'position',[0 0 17 13])
fa_field = 'fga_tmn';
for cc = 1:nc
cellid = pop_res_on(cc).cellid;
assert(cellid == pop_res_ms(cc).cellid)

ax = pop_res_on(cc).dv_axis';
xx = linspace(ax(1),ax(end),100);

%fa_on = mean(pop_res_on(cc).(fa_field));
fa_on = mean(pop_res_on(cc).fgta_resid(pop_res_on(1).t0s>.5,:));
%fa_ms = pop_res_ms(cc).(fa_field);
ms_goodtind = pop_res_on(cc).t0s>.25
fa_ms = mean(pop_res_ms(cc).fgta_resid(ms_goodtind,:));

subplot(ncols, nrows, cc);
plot(ax,fa_on,'.-k')
hold on
plot(ax,fa_ms,'.-r')
title(cellid)
%ylim([0 1])

drawnow
%pause(.5)
end
%%
mismatches = sign(slope_on ) ~= sign(slope_ms);
bad_ix = find(mismatches);
ax = [subplot(221) subplot(222) subplot(223) subplot(224)];
for ii = 1:length(bad_ix)
    ix = bad_ix(ii);
%ix = 84
fh = figure; clf
set(fh, 'position', [0 20 8 4])


fgta_line_plot(pop_res_on(ix), 'ax', subplot(221))
xlabel('time from onset')
ylabel('FR (z-score)')

fgta_line_plot(pop_res_ms(ix), 'ax', subplot(222))
xlabel('time from switch')
ylabel('FR (z-score)')

fgta_plot_tuning(pop_res_on(ix), 'ax', subplot(223))
ylabel('FR (z-score)')

fgta_plot_tuning(pop_res_ms(ix), 'ax', subplot(224))
ylabel('FR (z-score)')
drawnow
pause(.5)
end
%%
xval_id = 18362;
res_xv1 = psg.res_fn(xval_id,[],[],0)