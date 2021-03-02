alignment   = 'stimstart'; %'cpokeout'
lag         = 0.0;            %0.2;
t0s         = .1:0.025:1-lag;n_dv_bins   = 100;          % triggers rebinning, but not recompiling 
krn_width   = 0.1;
force_frdv  = 1;            % keep as one if rebinning
force_bin   = 0;
force_dv    = 1;
norm_type   = 'none';
ops.mask_after_stim_firing    = 1;
ops.mask_stim_firing          = 0;
ops.average_a_vals            = 1;
ops.average_fr                = 1;
ops.min_fr_mod                = 1;
ops.fit_time_sigmoid          = 0;
ops.use_nresidual             = 1;
ops.plot_res                  = 3;

excellid = 18181;

frbins      =  0:.5:100;
dp          = set_dyn_path;
direction   = 'backward';
frbins      =  0:1:50;
krn_type    = 'halfgauss';
%% 
thefit = fit_rat_analytical('H037','results_dir',dp.model_fits_dir)
%%
data = dyn_cell_packager(excellid);
[~, vec_data, ~,~] = get_behavior_data(dp.spikes_dir, ...
    data.cellid, data.sessid,ops);
[model, constant_x, allsame] = get_data_model_p(data, vec_data);
[align_strs, align_args] = dyn_align_LUT;
align_ind = strmatch(alignment,align_strs,'exact');
nanfrates   = isnan(data.frate{align_ind});
%%

%% re-analyze cell 18181 
dt = .1;
max_t = 1;
min_t = dt;
t0s         = min_t:dt:max_t-lag;  % triggers rebinning and recompiling


n_dv_bins = 20;
%n_dv_bins = [-10 -6:.1:6 10];
end_mask_s = 0.1;
which_switch = 'model';
res = dyn_fr_dv_map(excellid, 't0s', t0s, 'model', model, 'lag', lag, ...
    'frbins', frbins, 'alignment', alignment,...
    'krn_width', krn_width, 'krn_type', krn_type,...
    'norm_type', norm_type,'shuffle_trials',0,...
    'n_dv_bins',n_dv_bins,'end_mask_s',end_mask_s);

mvmnsz = 1;
plot_cell_map(res,'fh',1,'mvmnsz',mvmnsz)

%%
figure(10); clf
pa = squeeze(sum(res.pjoint_tfa,2));
subplot(211)
imagesc(pa,'x',res.dv_axis)
subplot(212); cla
plot(res.dv_axis,pa(end,:))
hold on
plot(res.dv_axis,mean(pa))
%%
res_shuff = dyn_fr_dv_map(excellid, 't0s', t0s, 'model', model, 'lag', lag, ...
    'frbins', frbins, 'alignment', alignment,...
    'krn_width', krn_width, 'krn_type', krn_type,...
    'norm_type', norm_type,'shuffle_trials',1,...
    'n_dv_bins',n_dv_bins,'end_mask_s',end_mask_s);

plot_cell_map(res_shuff,'fh',2,'mvmnsz',mvmnsz)

%% switch analysis
switch_t0s = [-.45:.025:.75]
switch_res = dyn_fr_dv_map(excellid, 't0s', switch_t0s, ...
    'model', model, 'lag', lag, ...
    'frbins', frbins, 'alignment', alignment,...
    'krn_width', krn_width, 'krn_type', krn_type,...
    'norm_type', norm_type,'shuffle_trials',0,...
    'n_dv_bins',[-8:.5:8],'end_mask_s',0,...
    'which_switch',which_switch);
plot_cell_map(switch_res,'fh',10)
%%
[x, frbins, Pjoints, fr_given_as, fr_var_given_as, a_given_frs] = ...
    fr_dv_switch(excellid, sw_t0s, ops,'lag', lag, ...
    'frbins', frbins, 'alignment', alignment,...
    'krn_width', krn_width, 'krn_type', krn_type,...
    'norm_type', norm_type, 'which_switch',which_switch);
%% analyze shuffle data for same cell

shuff_res = dyn_fr_dv_map(excellid, 't0s', switch_t0s, 'model', model, 'lag', lag, ...
    'frbins', frbins, 'alignment', alignment,...
    'krn_width', krn_width, 'krn_type', krn_type,...
    'norm_type', norm_type,'shuffle_trials',1,...
    'n_dv_bins',n_dv_bins,'end_mask_s',end_mask_s,...
    'which_switch',which_switch);

plot_cell_map(shuff_res,'fh',2)
%% analyze shuffle data for same cell
  % triggers rebinning and recompiling

fakefrates  = nan(size(data.frate{align_ind}));
ft          = data.frate_t{align_ind};
ft_         = [ft ft(end)+diff(ft(1:2))];
b = .025;

b= .01;

mu = 0;
ops.reload    = 0;
ops.ratname   = data.ratname;

b = 0;
b = 1;
fk = 2;

n = 10;
lc = 'k';

fnum = 5;

dv_bins = [-2.5:.25:2.5 ]
%dv_bins = 20;
bs = [.1 .25 .5 .75 1 2 10 100];
%bs = [.1 .25 1 ];
figure(fnum); clf
%%
bs = 100;
for b = bs
use_switch = 0;
switch 9
    case 10
        minf = -1.6;
        rangef = 1.6-minf;
        faketuning = @(dv,ft) minf + rangef*glmval([0 b]',dv(:),'logit');
        n = 0;
        fprintf('using the sign of the accumulator')
        lc = find(bs==b)/length(bs)*[1 0 1];
        if b > 1
            lc = 'r';
        end
    case 9
        faketuning = @(dv,ft) b*sign(dv(:));
        n = 0;
        tstr = 'sign(a)'
        fprintf('using the sign of the accumulator')
        lc = 'c';
        
    case 5 % the accumulator
        b = 1;
        faketuning = @(dv,ft) b*dv(:);
        n = 0;
        if b== 1
            tstr = 'a';
        else
            tstr = fprintf('%i a', tstr)
        end
        lc = 'k';
    case 0 % a tuning, no t mod
        faketuning = @(dv,ft) fk*ft(:) + ...
            20.*glmval([0 b]',[(dv(:)-b.*(ft(:)>0))],'logit') ;
    case 1 % a tuning, negative t mod
        faketuning = @(dv,ft) fk*ft(:) + ...
            20.*glmval([0 b]',[(dv(:)-b.*(ft(:)>0)).*(2-ft(:)).*.25],'logit') ;
    case 2 % a tuning, positive t mod
        faketuning = @(dv,ft) fk*ft(:) + ...
            20.*glmval([0 b]',[(dv(:)-b.*(ft(:)>0)).*(ft(:)).*.25],'logit') ;
    case 3 % no a tuning
        faketuning = @(dv,ft) fk*ft(:) ;
    case 4 % exponential tuning to accumulator
        k2= 1; a = 30;n=0;
        faketuning = @(dv,ft) a+k2.*exp(b*dv(:))+b.*ft(:);

  
    case 6 % steep accumulator + ramp
        beta = 5;
        faketuning = @(dv,ft) 1*ft(:)+10*glmval([0 beta]',[dv(:)],...
            'logit');
    case 7 % only choice tuned
        faketuning = @(dv,ft) 5+ft(:)+10*ft(:).^.25  * ...
            sign(dv(find(~isnan(dv),1,'last')));
    case 8 % choice tuned + a tuning
        faketuning = @(dv,ft) 5+ft(:)+5*ft(:).^.25  * ...
            sign(dv(find(~isnan(dv),1,'last'))) + ...
            1.*glmval([0 b]',[(dv(:)-b.*(ft(:)>0))],'logit');
    

end

fakenoisyfr = @(dv,ft) faketuning(dv,ft) + ...
    movmean(ft(:).*n.*randn(size(ft(:))),10);

for tt = 1:length(data.trials.gamma)
    
    this_lpulses = histcounts(data.trials.lpulses{tt},ft_);
    this_rpulses = histcounts(data.trials.rpulses{tt},ft_);
    pulse_diff   = cumsum(this_rpulses) - cumsum(this_lpulses);
    model_mean   = model(tt).posterior.mean;
    mm_interp    = interp1(model(tt).posterior.T, model_mean, ft);
    inputmean    = mm_interp;
    for ff = 1:length(ft)-1
        if ft(ff) <= 0 
            inputmean(ff) = 0;
        else
            this_ind = model(tt).posterior.T > ft(ff) & ...
                model(tt).posterior.T <= ft(ff+1);
            inputmean(ff) = mean(model_mean(this_ind));
        end
    end
    inputmean(ft<=0.0) = 0;
    %fakefrates(tt,:) = movmean(a+k*exprnd(pulse_diff)+b.*ft,20);
    fakefrates(tt,:) = fakenoisyfr(inputmean,ft);
end
fakefrates(nanfrates) = nan;



figure(110); clf
subplot(221)
imagesc(fakefrates)

subplot(223)
plot(ft, nanmean(fakefrates),...
    'color','k')
hold on
plot(ft, nanmean(fakefrates(data.trials.rat_dir==1,:)),...
    'color',dp.right_color,'linewidth',1.5)
plot(ft, nanmean(fakefrates(data.trials.rat_dir==-1,:)),...
    'color',dp.left_color,'linewidth',1.5)

title('PSTH')
xlabel('time (s)')
ylabel('firing rate')   
box off
xlim([-.5 2])


dfr = .5;
nfr = 50;
frbins = linspace(min(fakefrates(:)), max(fakefrates(:)),nfr);
switch use_switch
case 1
fakeres = dyn_fr_dv_map(excellid, 't0s', switch_t0s, ...
    'model', model, 'trialnums', 1:length(model),...
    'which_trials', 1:length(model),...
    'lag', lag, ...
    'frbins', frbins, 'alignment', alignment,...
    'krn_width', krn_width, 'krn_type', krn_type,...
    'norm_type', norm_type,'shuffle_trials',0,...
    'n_dv_bins',dv_bins,'end_mask_s',0,...
    'frates', fakefrates,...
    'which_switch',which_switch);
case 0
fakeres = dyn_fr_dv_map(excellid, 't0s', t0s, ...
    'model', model, 'trialnums', 1:length(model),...
    'which_trials', 1:length(model),...
    'lag', lag, ...
    'frbins', frbins, 'alignment', alignment,...
    'krn_width', krn_width, 'krn_type', krn_type,...
    'norm_type', norm_type,'shuffle_trials',0,...
    'n_dv_bins',dv_bins,'end_mask_s',0,...
    'frates', fakefrates);
end

figure(110); 
ax = subplot(224); cla
dv_axis = fakeres.dv_axis;
fakecurve = ones(size(dv_axis));
for ff = 1:length(dv_axis)
    fakecurve(ff) = faketuning(dv_axis(ff),1);
end
plot(ax,dv_axis,fakecurve,'-','linewidth',2)
ylim([min(fakecurve) max(fakecurve)]+[-.1 .1])

% end
% faketest = faketuning(dv_axis,ones(size(dv_axis)));
% plot(ax,dv_axis,faketuning(dv_axis,ones(size(dv_axis))),'linewidth',2)


ax = subplot(222); cla
dv_axis = fakeres.dv_axis;
cm = flipud(colormapLinear(dp.left_color,length(dv_axis)/2));
cm = [cm; colormapLinear(dp.right_color,length(dv_axis)/2)];
cm(length(dv_axis)/2,:) = [];
set(ax,'ColorOrder',cm,'NextPlot','ReplaceChildren')
hold(ax,'on')
for ff = 1:length(dv_axis)
    plot(ax,ft(:),faketuning(ones(size(ft))*dv_axis(ff),ft),...
        'linewidth',2)

end 
plot(ax,ft(:),faketuning(ones(size(ft))*dv_axis(1),ft),...
    'linewidth',2,'color',cm(1,:))
plot(ax,ft(:),faketuning(ones(size(ft))*dv_axis(end),ft),...
    'linewidth',2,'color',cm(end,:))


% [xx yy] = meshgrid(dv_axis,fakeres.t0s);
% tune_grid = reshape(faketuning(xx,yy),size(xx,1),size(xx,2));
% 
% plot(ax,fakeres.t0s,movmean(tune_grid,1),'linewidth',2)

%legend('t_{end}','t_0','location','southeast')
box off
xlabel('accumulation value (a)')
ylabel('firing rate')
title('generative tuning')

axh = plot_cell_map(fakeres,'fh',fnum,'linecolor',lc)
title(axh(end), tstr)
end
%%
mvres = 5;

ax = subplot(122); cla
dv_axis = fakeres.dv_axis;
cm = flipud(colormapLinear(dp.left_color,length(dv_axis)/2));
cm = [cm; colormapLinear(dp.right_color,length(dv_axis)/2)];

[xx yy] = meshgrid(dv_axis,fakeres.t0s);
tune_grid = reshape(faketuning(xx,yy),size(xx,1),size(xx,2));
set(ax,'ColorOrder',cm,'NextPlot','ReplaceChildren')
plot(ax,t0s,movmean(tune_grid,mvres),'linewidth',2)

%legend('t_{end}','t_0','location','southeast')
box off
xlabel('accumulation value (a)')
ylabel('firing rate')
title('generative tuning')

%%

fakeres.cellid = 99999;
plot_cell(fakeres)
set(figure(1),'position',[ 2.7778    3.3333   11.2083    7.3889])

%%  
fr_given_ta_cells = cat(3,res.fr_given_ta);



%% compute tuning curves for all cells
res = build_tuning_curves();

%% plot rank 1 decomp for all cells
%%
% cm = flipud(colormapLinear(dp.left_color,50));
% cm = [cm; colormapLinear(dp.right_color,50)];


figure(13); clf

pref_l = mean(res.tuning_fa(res.dv_axis>0,:)) < mean(res.tuning_fa(res.dv_axis<0,:));

good_cells = res.cell_num;
bad_cells  = find(~res.cell_dex);
%which_cells = which_cells(1:200)
ngood = length(good_cells);
do_norm = 1
if ~do_norm
    plot_fa = res.tuning_fa./res.tuning_alpha;
    plot_fr = res.tuning_fr./res.tuning_beta;
else
    plot_fa = res.tuning_fa(:,:);
    plot_fr = res.tuning_fr(:,:);
end

plot_fa(:,pref_l) = flipud(plot_fa(:,pref_l));

plot_fafr = plot_fa .* mean(plot_fr);

mtrange = [percentile(plot_fr(:),[1 99])]
rarange = [percentile(plot_fa(:),[1 99])]

sort_cells = [good_cells; bad_cells];
[~, sort_cells] = sort(res.fr_modulation,'descend');


good_color = [.8 .9 .8];
bad_color = [.9 .8 .9];
subplot(231);
imagesc(plot_fr(:,sort_cells)','x',res.t0s)
hold on;
plot(xlim,ngood.*[1 1],'m')
ylabel('cell #')
xlabel('time from stim on (s)')
if do_norm
    title('m(t)=u1*s1*range(v1)')
else
    title('m(t)=u1')
end
caxis(mtrange)
cb = colorbar('southoutside')
subplot(234);
plot(res.t0s,plot_fr(:,bad_cells)','color',bad_color)
hold on
plot(res.t0s,plot_fr(:,good_cells)','color',good_color)

plot(res.t0s,nanmean(plot_fr(:,bad_cells)'),'color',bad_color.^10,'linewidth',2)

plot(res.t0s,nanmean(plot_fr(:,good_cells)'),'color',good_color.^10,'linewidth',2)

plot(xlim,[0 0],'k--')
ylim(mtrange)
ylabel('f.r. modulation (Hz)')
xlabel('time from stim on (s)')
title('m(t)')

subplot(232);
imagesc((plot_fa(:,sort_cells)'),'x',res.dv_axis)
hold on;
plot(xlim,ngood.*[1 1],'m')
ylabel('cell #')
xlabel('accumulation value (a)')
caxis(rarange)
if do_norm
    title('r(a)=v1/range(v1)')
else
    title('r(a)=v1')
end
cb = colorbar('southoutside')


subplot(233)
imagesc((plot_fafr(:,sort_cells)'),'x',res.dv_axis)
hold on;
plot(xlim,ngood.*[1 1],'m')
title('<m(t)> * r(a)')
cb = colorbar('southoutside')
colormap(cm)

subplot(235);
plot(res.dv_axis,plot_fa(:,bad_cells)','color',bad_color)
hold on
plot(res.dv_axis,plot_fa(:,good_cells)','color',good_color)
plot(res.dv_axis,mean(plot_fa(:,good_cells)'),'color',good_color.^10,'linewidth',2)
plot(res.dv_axis,mean(plot_fa(:,bad_cells)'),'color',bad_color.^10,'linewidth',2)
ylabel('normalized FR')
xlabel('accumulation value (a)')
title('r(a)')
ylim(rarange)


subplot(236)

plot(res.dv_axis,plot_fafr(:,good_cells)','color',good_color)
hold on
plot(res.dv_axis,plot_fafr(:,bad_cells)','color',bad_color)
plot(res.dv_axis,nanmean(plot_fafr(:,good_cells)'),'color',good_color.^10,'linewidth',2)
plot(res.dv_axis,nanmean(plot_fafr(:,bad_cells)'),'color',bad_color.^10,'linewidth',2)
title('<m(t)> * r(a)')
text(4,min(ylim)+.1*range(ylim),'good','color',good_color.^10,'fontsize',18)
text(4,min(ylim)+.2*range(ylim),'bad','color',bad_color.^10,'fontsize',16)

colormap(subplot(231),colormapRedBlue)
caxis(subplot(231),[-1 1].*10)
colormap(subplot(232),colormapRedBlue)
caxis(subplot(232),[-1 1].*1)
colormap(subplot(233),colormapRedBlue)
caxis(subplot(233),[-1 1].*1)

%%
pref_l = mean(res.tuning_fa(res.dv_axis>0,:)) < ...
    mean(res.tuning_fa(res.dv_axis<0,:));

figure(2); clf
which_cells = find(res.cell_dex);
plot_val = res.fga_nta_cell;
%plot_val = squeeze(nanmean(res.fga_cell_residual,1));
%plot_val = squeeze(nanmean(res.fga_cell_nresidual,1));
plot_val(:,pref_l) = flipud(plot_val(:,pref_l));

subplot(121)
x = res.dv_axis;
this_plot = plot_val(:,which_cells);
plot(x,this_plot,'color',[1 1 1].*.5)
hold on
plot(x,nanmean(this_plot'),'k','linewidth',2)
which_cells = find(~res.cell_dex);
subplot(122)
x = res.dv_axis;
this_plot = plot_val(:,which_cells);
plot(x,this_plot,'color',[1 1 1].*.5)
hold on
plot(x,nanmean(this_plot'),'k','linewidth',2)
%%
figure(100); clf
histogram(res.rank1_variance(res.cell_dex),'facecolor','k')
hold on
histogram(res.rank1_variance(~res.cell_dex),'facecolor','m')

nanmean(res.rank1_variance(res.cell_dex))
nanmean(res.rank1_variance(~res.cell_dex))
%%
cellid = 18181;

%{ 
there's one problem with shuffling, which is that not all the stimuli line
up with the firing rates. 
I could use unmasked trials, but the firing rate
statistics are different after cpoke out.
Instead, the smarter thing to do might be shuffle by quintile or something
like that.
%}

shuff_exres = dyn_fr_dv_map(cellid, t0s, n_dv_bins, ops,...
    'lag', lag, 'krn_width', krn_width, 'alignment', alignment, ...
    'var_weight', false, 'force_frdv',1,'force_bin',1,...
    'force_dv', 1,'norm_type',norm_type,'shuffle_trials',1,...
    'save_map',0);
plot_cell(shuff_exres)

exres = dyn_fr_dv_map(cellid, t0s, n_dv_bins, ops,...
    'lag', lag, 'krn_width', krn_width, 'alignment', alignment, ...
    'var_weight', false, 'force_frdv',1,'force_bin',1,...
    'force_dv', 1,'norm_type',norm_type);
plot_cell(exres)
%%
cellid = 18181;

%%


%%
Tsort = sort(data.trials.T);
maxT = Tsort(300);
maxT = 1.25
goodt = fakeres.t0s < maxT;
[u,s,v] = svd(fakeres.fga_cell_residual(goodt,:));
s2 = s(1:rank,1:rank);
%s2(2:end) = 0;
tuning_cell = u(:,1:rank)*s2*v(:,1:rank)';
alpha   = 1/range(v(:,1));
beta    = s2(1,1)/alpha;
tuning_fr = u(:,1)*beta;
tuning_fa = v(:,1)*alpha;
s_squared = diag(s).^2;
rank1_variance = s_squared(1)./sum(s_squared);
tuning_alpha = alpha;
tuning_beta = beta;

figure(153); clf
subplot(231)
plot(fakeres.t0s(goodt),tuning_fr)
subplot(232)
plot(fakeres.dv_axis,tuning_fa)
subplot(233)
imagesc(tuning_cell','x',fakeres.t0s(goodt),'y',fakeres.dv_axis)
subplot(234)
imagesc(fakeres.fga_residual(goodt,:)','x',fakeres.t0s(goodt),'y',fakeres.dv_axis)

%%
figure(150); clf; 
[f, x] = ecdf(data.trials.T,'function','survivor');
stairs([0; x],[1; f].*length(data.trials.T))



%% re-analyze cell 18181 
excellid = 18181;
shuff_trials = 0;
[x, frbins, Pjoints, fr_given_as, fr_var_given_as, a_given_frs] = ...
    dyn_compile_dv2(excellid, t0s, ops,'lag', lag, ...
    'frbins', frbins, 'alignment', alignment,...
    'krn_width', krn_width, 'krn_type', krn_type,...
    'norm_type', norm_type,'shuffle_trials',shuff_trials);

%%
fr_given_as = fakeres.fga_cell;
fga_cell_residual = fr_given_as - nanmean(fr_given_as,2);

cm = flipud(colormapLinear(dp.left_color,50));
cm = [cm; colormapLinear(dp.right_color,50)];

mvres = 1;
figure(11); clf
ax1 = subplot(221)
cm = color_set(size(fga_cell_residual(1:mvres:end,:),2));
set(ax1,'ColorOrder',cm,'NextPlot','ReplaceChildren')
plot(ax1,t0s,movmean(fga_cell_residual,mvres))
hold(ax1,'on')
%plot(x,nanmean(fr_given_as),'color','m','linewidth',2)
cb = colorbar('northoutside')
title(cb,'accumulation value (a)')
colormap(ax1,cm)
xlabel('time (s)')
ylabel('firing rate')
axis('tight')

ax2 = subplot(222)
cm = colormapLinear([ 1 1 1].*0, size(fga_cell_residual,1)).^.8;
set(ax2,'ColorOrder',cm,'NextPlot','ReplaceChildren')
plot(ax2,x,movmean(fga_cell_residual,mvres))
hold(ax2,'on')
plot(x,nanmean(fga_cell_residual),'color','m','linewidth',2)
cb = colorbar('northoutside')
title(cb,'time (s)')
colormap(ax2,cm)
xlabel('accumulation value (a)')
ylabel('firing rate')
axis('tight')


linkaxes([ax1, ax2],'y')

if shuff_trials
    suptitle(['cell ' num2str(excellid) ' SHUFFLE'])
else
    suptitle(['cell ' num2str(excellid) ' REAL'])
end
rank = 1;
[u,s,v] = svd(fga_cell_residual);
s2 = s(1:rank,1:rank);
%s2(2:end) = 0;
tuning_cell = u(:,1:rank)*s2*v(:,1:rank)';
alpha   = 1/range(v(:,1));
beta    = s2(1,1)/alpha;
tuning_fr = u(:,1)*beta;
tuning_fa = v(:,1)*alpha;
s_squared = diag(s).^2;
rank1_variance = s_squared(1)./sum(s_squared);
tuning_alpha = alpha;
tuning_beta = beta;
% sign of
if median(tuning_fr) < 0
    tuning_fr = -u(:,1)*beta;
    tuning_fa = -v(:,1)*alpha;
end
subplot(223);
plot(t0s,tuning_fr,'k','linewidth',1.5)
hold on
plot(xlim,[ 0 0],'k')
box off
ylabel('m(t)')
xlabel('time (s)')
subplot(224);
plot(x,tuning_fa,'k','linewidth',1.5)
hold on
plot(xlim,[ 0 0],'k')

box off
ylabel('r(a)')
xlabel('accumulation value (a)')

%cellid = [18181 18185 17799 17784 16875 18172 17870 16847 17855 17803 18466]
%cellid = [-1 -2 -3 -4 -5];
%cellids = [-1 -11 -2 -20 -21 -22 -23 -24 -3 -30 -31 -32];
%%
excellid = 16857;
exres = dyn_fr_dv_map(excellid, t0s, n_dv_bins, ops, ...
    'lag', lag, 'krn_width', krn_width, 'alignment', alignment, ...
    'var_weight', false, 'force_frdv',force_frdv,'force_bin',force_bin, ...
    'force_dv', force_dv,'norm_type',norm_type,'frbins',frbins);
 
plot_cell(exres)
%%
figure(13); clf
subplot(221);
imagesc(exres.tuning_fr','x',exres.t0s)
colormap(bone)
subplot(223);
plot(exres.t0s,exres.tuning_fr','color',[.5 .5 .5])
hold on
plot(exres.t0s,mean(exres.tuning_fr'),'color',[.5 .5 .5].*0,'linewidth',2)
subplot(222);
imagesc(exres.tuning_fa','x',exres.dv_axis)
subplot(224);
plot(exres.dv_axis,exres.tuning_fa','color',[.5 .5 .5])
hold on
plot(exres.dv_axis,mean(exres.tuning_fa'),'color',[.5 .5 .5].*0,'linewidth',2)



%%
res = build_tuning_curves('H037');
%%
ii = 5
which_cell = find(res.cell_dex)
plot_cell(res,which_cell(ii))
%%
res2 = build_tuning_curves('H066');
%%
which_cell = find(res2.cellid==18181);
plot_cell(res2,which_cell)
%%
res3 = build_tuning_curves('H191');


%%
figure(12); clf
subplot(121)
plot(exres.tuning_fr)
subplot(122)
plot(exres.tuning_fa)

%%
excellid = [18181 16857 17784];
exres = dyn_fr_dv_map(excellid, t0s, n_dv_bins, ops,...
    'lag', lag, 'krn_width', krn_width, 'alignment', alignment, ...
    'var_weight', false, 'force_frdv',force_frdv,'force_bin',force_bin,...
    'force_dv', force_dv,'norm_type',norm_type,'frbins',frbins);
%% 
plot_population(exres)
%%

[x, frbins, Pjoints, fr_given_as, fr_var_given_as, a_given_frs] = ...
    dyn_compile_binned_database(excellid, t0s, n_dv_bins, ops,...
    'lag', lag, 'krn_width', krn_width, 'fr_dt', [], 'alignment', alignment,...
    'direction', direction, 'frbins', frbins, 'krn_type', krn_type, 'norm_type', norm_type,...
    'n_iter',1, 'param_scale_num', 1, ...
    'param_scale_factor', 1, 'force_bin', force_bin,'force_dv',force_dv,...
    'datadir',dp.celldat_dir);%, 'bootstrap', bootstrap);

%%
excellid = 18181;
which_switch = 'model';
sw_t0s = [-.45:.025:.75]
[x, frbins, Pjoints, fr_given_as, fr_var_given_as, a_given_frs] = ...
    fr_dv_switch(excellid, sw_t0s, ops,'lag', lag, ...
    'frbins', frbins, 'alignment', alignment,...
    'krn_width', krn_width, 'krn_type', krn_type,...
    'norm_type', norm_type, 'which_switch',which_switch);
%%
figure(13); clf
plot(x,fr_given_as,'color',[1 1 1].*.7)
hold on
plot(x,nanmean(fr_given_as(sw_t0s<0,:)),'color','k','linewidth',2)
plot(x,nanmean(fr_given_as(sw_t0s>0,:)),'color','k','linewidth',2)
%%
cm = flipud(colormapLinear(dp.left_color,size(fr_given_as,2)/2));
cm = [cm; colormapLinear(dp.right_color,size(fr_given_as,2)/2)];
figure(14); clf
ax = axes
%cm =color_set(size(fr_given_as,2));
set(ax,'ColorOrder',cm,'NextPlot','ReplaceChildren')
plot(ax,sw_t0s,fr_given_as)
xlabel(ax,['time from ' which_switch ' switch'])
colormap(cm);
cb = colorbar('northoutside')
title(cb,'accumulated value (a)')
%%
[x, frbins, Pjoints, fr_given_as, fr_var_given_as, a_given_frs] = ...
    dyn_compile_dv2(excellid, t0s, ops,'lag', lag, ...
    'frbins', frbins, 'alignment', alignment,...
    'krn_width', krn_width, 'krn_type', krn_type,...
    'norm_type', norm_type);
figure(11); clf
plot(x,fr_given_as,'color',[1 1 1].*.7)
hold on
plot(x,nanmean(fr_given_as),'color','k','linewidth',2)
%%
figure; imagesc(a_given_frs,'x',frbins,'y',t0s)


%%
switch_t0s = [-.45:.025:.75]

