alignment   = 'stimstart'; %'cpokeout'
lag         = 0.0;            %0.2;
dt = .025;
t0s         = .1:dt:1-lag;n_dv_bins   = 100;          % triggers rebinning, but not recompiling
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
switch_t0s = [-.45:.025:.55];

p.compute_full  = 1;
p.compute_back  = 1;
p.compute_delta = 5;
p.save_new_ref  = 0;

%p.dt = .01;
%p.check_final   = 1;

%%
max_T = 2;
p.dt            = 1e-2  ;
t = 0:p.dt:max_T-p.dt;
p.n             = 100;
p.b             = 50;
p.da            = 1;
p.avals         = -p.b:p.da:p.b;
p.d_dex         = 6; % abs value of delta-mode to compute
p.error_tolerance = 0.0001;
p.a_sign = 1
%%
dv_bins     = [-5 : .25 : 5];
frbins      =  0:.5:100;
dp          = set_dyn_path;
direction   = 'backward';
%frbins      =  0:1:50;
krn_type    = 'halfgauss';
%%
figure(1); clf
xx = -10:.1:10;
yy = glmval([0 .01]',xx(:),'logit');
yyn = yy - min(yy) ;
yyn = yyn ./ max(yyn);
[b4] = fit_four_param_psycho(xx(:), yyn(:), 0)
b = glmfit(xx(:),yyn(:),'binomial')
subplot(211)
plot(xx(:),yyn(:),'o')
hold on
plot(xx(:),fourParamPsychometric(b4,xx(:)),'m','linewidth',1);
subplot(212)
yhat = glmval(b,xx(:),'logit');
plot(xx(:),yyn(:),'.k')
hold on
plot(xx,yhat','m')
%%
cellid      = 18181;
celldata    = dyn_cell_packager(cellid)
%%
ratname     = celldata.ratname;
sessid      = celldata.sessid;
%%
thefit = fit_rat_analytical(ratname,'data_dir', dp.behav_data_dir, ...
    'results_dir',dp.model_fits_dir);
%%
params_real = thefit.final;
param_names = {'\lambda', '\sigma_a', '\sigma_s', '\sigma_i', '\phi',  '\tau', 'bias', 'lapse'};
params_nonoise = params_real;%.* [1 0 0 0 1 1 0 0];
params_lessnoise = params_real;
params_nonoise([2 3 4]) = [0 0 0];
params_lessnoise([2 3 4]) = [eps .5 eps];

%%
[array_data, vec_data, ~,~] = get_behavior_data(dp.spikes_dir, ...
    celldata.cellid, celldata.sessid,ops);
[align_strs, align_args] = dyn_align_LUT;
alignment = 'stimstart' ;
align_ind = strmatch(alignment,align_strs,'exact');
nanfrates   = isnan(celldata.frate{align_ind});
cell_ft      = celldata.frate_t{align_ind};
%%

NT                  = length(celldata.trials.T);
asample_lessnoise         = nan(NT,numel(t));
pokedR_lessnoise     = nan(NT,1);
asample_nonoise         = nan(NT,numel(t));
pokedR_nonoise     = nan(NT,1);
asample_ratnoise         = nan(NT,numel(t));
pokedR_ratnoise     = nan(NT,1);
match_choice    = 0;
tnums = celldata.trials.trialnums;
%%
sample_params = params_lessnoise;
for  tn = 1:NT;
    %%
    if mod(tn,10)==0
        disp(tn)
    end
    trial.T          = celldata.trials.T(tn);
    trial.leftbups   = celldata.trials.lpulses{tn};
    trial.rightbups  = celldata.trials.rpulses{tn};
    trial.rat_dir    = celldata.trials.rat_dir(tn);
    
    particle_lessnoise = compute_particles(trial, params_lessnoise, p, []);
    particle_nonoise = compute_particles(trial, params_nonoise, p, []);
    particle_ratnoise = compute_particles(trial, params_real, p, []);
    
    lessnoise_sample   = particle_lessnoise.a(1,:);
    nonoise_sample     = particle_nonoise.a(1,:);
    thistind = 1:length(lessnoise_sample);
    tvec = (0:length(thistind)-1)/p.dt;
    asample_lessnoise(tn,thistind)   = lessnoise_sample;
    asample_nonoise(tn,thistind)   = nonoise_sample;
    pokedR_lessnoise(tn) = lessnoise_sample(end) > 0;
    pokedR_nonoise(tn) = nonoise_sample(end) > 0;
    

    pchoice         = particle_ratnoise.a(:,end) > 0;
    samechoice      = find(pchoice == pokedR_lessnoise(tn));
    samechoiceind   = samechoice(randperm(length(samechoice)));
    ratnoise_sample   = particle_lessnoise.a(samechoiceind(1),:);
    asample_ratnoise(tn,thistind)   = ratnoise_sample;
    pokedR_ratnoise(tn) = ratnoise_sample(end) > 0;
    
end
%%
model_real   = get_data_model_p( sessid,tnums)

%%
avals = model_real(1).posterior.avals;
model_nonoise = model_real;
avals_edges = avals(1:end-1) + diff(avals)/2;


for tn = 1:length(model_nonoise)
    M = model_nonoise(tn).posterior.pdf;
    [i j] = max(M,[],2);
    M = zeros(size(M));
%     [~, ~, j] = histcounts(asample(mm,~isnan(asample(mm,:))),avals_edges);
%     j = j+1;
%     for ii = 1:length(j)
%         M(ii,j(ii)) = 10;
%     end
    var = .1;
    this_sample = asample_nonoise(tn,:);
    M = nan(numel(this_sample), numel(avals));
    for ii = 1:size(M,1)
        M(ii,:) = normpdf(avals,asample_nonoise(tn,ii),var);
        M(ii, :) = M(ii,:) ;
    end
    model_nonoise(tn).posterior.pdf = M;
    MT = model_nonoise(tn).posterior.T;
    model_nonoise(tn).posterior.T = t(1:numel(this_sample));
end

%%
tnums = celldata.trials.trialnums;
model_ratnoise   = compute_model_mean(ratname, sessid, ...
    'fake_params', params_real, 'fake_choices', pokedR_ratnoise,...
    'which_trials',tnums,'dt',p.dt)
%%
savename = fullfile(dp.model_dir,'simulate_model_ratnoise');
save(savename, 'asample_ratnoise','model_ratnoise','ratname','sessid','t');

%% TEMPORARY
fr_gen = sign(asample_ratnoise);
ft = t;
which_switch = 'model';
fn = 3;
dv_bins = [-4:2:4]
switch_t0s = -.55:.025:.55;
lag = 0;
switch fn
    case 1
        hits_only = 0;
        errs_only = 0;
        
    case 2
        hits_only = 1;
        errs_only = 0;
        %fr_gen = fr_gen(celldata.trials.hit==1,:);
    case 3
        hits_only = 0;
        errs_only = 1;
        %fr_gen = fr_gen(~celldata.trials.hit==0,:);
end

res_switch = dyn_fr_dv_map(18181, 't0s', switch_t0s, ...
    'frates', fr_gen, 'ft', ft,...
    'model', model_ratnoise, 'data', celldata,...
    'lag', lag, ...
    'frbins', frbins, 'alignment', alignment,...
    'krn_width', krn_width, 'krn_type', krn_type,...
    'norm_type', norm_type,'shuffle_trials',0,...
    'n_dv_bins',dv_bins,'end_mask_s',0,...
    'which_switch',which_switch,...
    'hits_only',hits_only,'errs_only',errs_only,...
    'do_svd',0);

%%
res_switch = compute_rank1_fgta_approx(res_switch);

figure(fn); clf
ax = axes;
fgta_line_plot(res_switch,'ax',ax)

figure(fn+1); clf
plot(res_switch.dv_axis,  res_switch.rank1_ra_n, '.-')

%%
lag = .1;
res_switch = dyn_fr_dv_map(18181, 't0s', switch_t0s, ...
    'model', model_real, 'data', celldata,...
    'lag', lag, ...
    'frbins', frbins, 'alignment', alignment,...
    'krn_width', krn_width, 'krn_type', krn_type,...
    'norm_type', norm_type,'shuffle_trials',0,...
    'n_dv_bins',dv_bins,'end_mask_s',0,...
    'which_switch',which_switch,...
    'hits_only',1,'errs_only',0);
%%
figure(100); clf

fgta_line_plot(res_switch,'linewidth',2)

%%
model_lessnoise   = compute_model_mean(ratname, sessid, ...
    'fake_params', params_lessnoise, 'fake_choices', pokedR_lessnoise,...
    'which_trials',tnums,'dt',p.dt)


%%
labels  = {'nonoise', 'lessnoise', 'ratnoise'};
dvlims = [-3 3]
figure(3); clf

tn = 100;
s(1) = subplot(311)
imagesc(model_nonoise(tn).posterior.pdf', ...
    'y', model_nonoise(tn).posterior.avals,...
    'x', model_nonoise(tn).posterior.T)
colormap(flipud(gray))
hold on
plot(t,asample_nonoise(tn,:),'r')
xlim([0 2.1])
caxis([0 1])
title(labels{1})
s(2) = subplot(312);
imagesc(model_lessnoise(tn).posterior.pdf', ...
    'y', model_lessnoise(tn).posterior.avals,...
    'x', model_lessnoise(tn).posterior.T)
hold on
plot(t,asample_lessnoise(tn,:),'r')
xlim([0 2.1])
caxis([0 1])
title(labels{2})
s(3) = subplot(313)
imagesc(model_ratnoise(tn).posterior.pdf', ...
    'y', model_ratnoise(tn).posterior.avals,...
    'x', model_ratnoise(tn).posterior.T)
hold on
plot(t,asample_ratnoise(tn,:),'r')
xlim([0 2.1])
linkaxes(s)
ylim(dvlims)
caxis([0 1])
ylabel('a')
xlabel('t')
title(labels{3})
%%
figure(1); clf
figure(2); clf

clrs    = [[0 0 0];  [.25 .75 .25]; [.75 .25 .75]; [.75 .75 .25]]; 
faketune = @(x, slope) 5+10*glmval([0 slope]', x, 'logit');
slope_in = [.1 1 10 100];
remove_noise = 0;
add_spike_noise = 1;

for ll = 1:length(labels)
    slope_out = nan(size(slope_in));
for ss = 1:length(slope_in)
    
asample     = eval(['asample_' labels{ll}]);
model       = eval(['model_' labels{ll}]);
fr_gen      = faketune(asample(:), slope_in(ss));
fr_gen      = reshape(fr_gen, size(asample,1), size(asample,2));
fr_in       = fr_gen;
ft          = t;
%%
if add_spike_noise
    spikes = fr_gen * p.dt > rand(size(fr_gen)) ;
    krndt   = .1;
    dt      = p.dt;
    krn     = my_gauss_krn(krndt,'causal',dt,1);
    x       = 0:krn.bin_size:max_T-krn.bin_size;
    fr_hat = nan(NT, numel(x));
    for tn = 1:NT
        [y ~] = spike_filter(0,t(spikes(tn,:)), krn.krn, ...
            'pre',t(1),'post',max_T,'kernel_bin_size', krn.bin_size);
        fr_hat(tn,:) = y;
        ind = find(nanfrates(tn,:),1);
        if ~isempty(ind)
            ind = find(t>cell_ft(ind),1);
            
            fr_hat(tn,ind:end) = nan;
        end
    end
    fr_in   = fr_hat;
    ft      = x;
    
    figure(100); clf
    tn = 30
    subplot(3,4,1)
    %avals = -10:.1:10;
    plot(avals',faketune(avals, slope_in(ss)),'r','linewidth',2)
    title('tuning')
    ylabel('fr')
    xlabel('a')
    %ylim([-.5 10.5])
    
    
    subplot(3,4,[2 3 4])
    plot(t,asample(tn,:),'linewidth',2)
    hold on
    plot(t,fr_gen(tn,:),'.m','linewidth',1,'markersize',12)
    
    plot(ft, fr_hat(tn,:),'-g','linewidth',1.5)
    this_spks = t(spikes(tn,:));
    if ~isempty(this_spks)
        plot(this_spks,11,'.k','markersize',10)
    end
    plot(xlim,[0 0],'-k')
    
    xlabel('time from stim on (s)')
    legend('a','real rate',  'est rate','spikes','location','eastoutside')
    title('example trial')
    
end
%%
if remove_noise
    model_ = model;
    for tn = 1:NT   
        M = model(tn).posterior.pdf;
        M_ = zeros(size(M));
        [i j] = max(M,[],2);
        for ii = 1:length(j)
            M_(ii,j(ii)) = 10;
        end
        model_(tn).posterior.pdf = M_;
    end
    %%
    %figure(10); clf
    figure(100);
    tn = 30
    subplot(3,1,2)
    imagesc(model_(tn).posterior.pdf','x',model_(tn).posterior.T,...
        'y',model_(tn).posterior.avals)
    colormap(flipud(bone))
    hold on
    plot(t,asample(tn,:),'r')
    ylim(dvlims)
    subplot(3,1,3)
    imagesc(model(tn).posterior.pdf','x',model(tn).posterior.T,...
        'y',model(tn).posterior.avals)
    colormap(flipud(bone))
    hold on
    plot(t,asample(tn,:),'r')
    model = model_;
    ylim(dvlims)
end
%%

%yy = glmval([0 .01]',xx(:),'logit');

frbins = linspace(nanmin(fr_in(:)),nanmax(fr_in(:)), 50);
res(ss) = dyn_fr_dv_map(cellid, 't0s', t0s, 'model', model, 'lag', lag, ...
        'frates', fr_in, 'ft', ft, ...
        'which_trials',1:length(model),'trialnums',1:length(model),...
        'frbins', frbins, 'alignment', alignment,...
        'krn_width', krn_width, 'krn_type', krn_type,...
        'norm_type', norm_type,'shuffle_trials',0,...
        'n_dv_bins',dv_bins,'end_mask_s',0,...
        'average_fr', 0, 'average_a_vals', 0);
    
%%
tune_field = 'fr_given_ta';
xx = res(ss).dv_axis;
dvlims = [-2.5 2.5];
fit_x = xx > dvlims(1) & xx < dvlims(2);
xx = xx(fit_x);
yy = mean(res(ss).(tune_field)(end-5:end,fit_x));
yyn = yy - min(yy) ;
yyn = yyn ./ max(yyn);
[b4] = fit_four_param_psycho(xx(:), yyn(:), 0)
slope_out(ss) = b4(3);

figure(1);
pnum = numel(slope_in)*(ss-1)+ll;
pnum = numel(labels)*(ss-1)+ll;
subplot(numel(slope_in),numel(labels),pnum)
plot(xx,yyn)
tuning = faketune(xx(:), slope_in(ss));
tuning = tuning - min(tuning(:));
tuning = tuning ./ max(tuning(:));

plot(xx(:), tuning, 'color',[.5 .5 .5].*0)
hold on
plot(xx(:),fourParamPsychometric(b4,xx(:)), 'color','m','linewidth',1);
plot(xx(:),yyn(:),'o','markeredgecolor',clrs(ll,:),...
    'markerfacecolor','w','markersize',7,'linewidth',2)

drawnow
xlim(dvlims+[-.1 .1])
ylim([-.1 1.1])
ylabel('norm fr')
xlabel('a')
title({labels{ll} [' slope = ' num2str(slope_in(ss))]})
%axis square
end
%%
figure(2); 
semilogx(slope_in, slope_out,'o-','linewidth',2,'color',clrs(ll,:))
hold on
end

xlabel('generative slope')
ylabel('recovered slope')
hl = legend(labels,'location','northwest')
%%
figure(4); 

%%
yyn = tuning;
[b4] = fit_four_param_psycho(xx(:), tuning(:), 0)


%% THIS IS WHERE WE CAN PLAY WITH THE VARIANCE
var = 2;
M = nan(numel(lessnoise_sample), numel(avals));
for ii = 1:size(M,1)
    M(ii,:) = normpdf(avals,lessnoise_sample(ii),var);
    M(ii, :) = M(ii,:) ; 
end
figure(1); clf
s(1) = subplot(311)
imagesc(M','x',t(thistind),'y',avals)
colormap(flipud(bone))
hold on
plot(t(thistind),lessnoise_sample,'r')
colorbar
caxis([0 .25])
s(2) = subplot(312)
imagesc(model_real(1).posterior.pdf')
colorbar
caxis([0 .25])
tt = 100;
subplot(313)
plot(avals,M(tt,:))
hold on
plot(avals,model_real(1).posterior.pdf(tt,:))
%%

%%
Mdvs = model_fake_noisy(1).posterior.avals;
MT = model_fake_noisy(1).posterior.T;
M = model_fake_noisy(1).posterior.pdf;
subplot(212)
imagesc(M','x',MT,'y',Mdvs)
hold on
ylim([-1 1].*3)
axis xy
plot(t,asample(1,:),'r')
%%

%%
pwr = 100;
model_fake = model_fake_noisy;
%%
for tn = 1 : NT
    %%
Mdvs    = model_fake_noisy(tn).posterior.avals;
M       = model_fake_noisy(tn).posterior.pdf;
MT      = model_fake_noisy(tn).posterior.T;
M       = M./sum(M,2);
M2      = M .^ pwr;
M2      = M2 ./ sum(M2,2) ;
oldmn   = sum(M.*Mdvs,2); 
newmn   = sum(M2.*Mdvs,2);
newvarx = (Mdvs - newmn).^2;
newvars = sum(M2.*newvarx,2);
newvar  = mean(newvars);

model_fake(tn).posterior.pdf = M2;

oldvarx = (Mdvs - newmn).^2;
oldvars = sum(M.*oldvarx,2);

figure(1); clf
subplot(211)
imagesc(M','x',MT,'y',Mdvs)
hold on
plot(MT,oldmn,'r','linewidth',1.5)
ylim([-1 1].*3)
axis xy
plot(t,asample(tn,:),'m')

%plot(MT,oldvars,'b')

subplot(212)
imagesc(M2','x',MT,'y',Mdvs)
axis xy
hold on
plot(MT,newmn,'r','linewidth',1.5)
plot(t,asample(tn,:),'c')
ylim([-1 1].*3)
colormap(flipud(bone))
drawnow
pause()
end
%%

    %%
plot_cell_map(res)
%%

figure(1); clf
subplot(211)
imagesc(M','x',MT,'y',Mdvs)
hold on
plot(MT,oldmn,'r','linewidth',1.5)
ylim([-1 1].*3)
axis xy

plot(MT,oldvars,'b')

subplot(212)
imagesc(M2','x',MT,'y',Mdvs)
axis xy
hold on
plot(MT,newmn,'r','linewidth',1.5)

ylim([-1 1].*3)
colormap(flipud(bone))


plot(MT,newvars*pwr,'b')
plot(MT,newvars,'g')

%%
avals = [-5:.25:5];


%%
%%
avals = [-5:.25:5];

avals = model_noisy(1).posterior.avals;
avals_edges = avals(1:end-1) + diff(avals)/2;


switch 2
    case 0
        params = params_real;
        model = model_real;
        match_choice = 1;
    case 1
        params = params_nonoise;
        model = model_nonoise;
        match_choice = 0;
        
    case 2
        model_nonoise = model_noisy;
        for tn = 1:length(model_nonoise)
            M = model_nonoise(tn).posterior.pdf;
            [i j] = max(M,[],2);
            M = zeros(size(M));
            [~, ~, j] = histcounts(asample(tn,~isnan(asample(tn,:))),avals_edges);
            j = j+1;
            for ii = 1:length(j)
                M(ii,j(ii)) = 10;
            end
            model_nonoise(tn).posterior.pdf = M;
        end
        params_nonoise([2 3 4]) = [eps eps eps];
        params = params_nonoise;
        model = model_nonoise;
end
%%
gen_fr = nan(NT,numel(t));

for tn = 1:NT;
    lessnoise_sample = asample(tn,:);
    gen_fr(tn,1:length(lessnoise_sample)) = faketuning(lessnoise_sample,tvec);
end

spikes = gen_fr * max_T * p.dt > rand(size(gen_fr)) ;
krndt   = .15;
dt = p.dt;
krn    = my_gauss_krn(krndt,'fullgauss',dt,1);
x = 0:krn.bin_size:max_T-krn.bin_size;
fr_hat = nan(NT, numel(x));
for tn = 1:NT
    [y x] = spike_filter(0,t(spikes(tn,:)), krn.krn, ...
        'pre',t(1),'post',max_T,'kernel_bin_size', krn.bin_size);
    fr_hat(tn,:) = y;
    ind = find(nanfrates(tn,:),1);
    if ~isempty(ind)
        ind = find(t>cell_ft(ind),1);
        
        fr_hat(tn,ind:end) = nan;
    end
end
%fr_hat(nanfrates) = nan
ft = t;
%%

tn = 200
%fr_hat = (asample);

figure(104); clf


subplot(221)
%avals = -10:.1:10;
plot(avals',faketuning(avals, ones(size(avals))),'r','linewidth',2)
title('tuning')
ylabel('fr')
xlabel('a')
%ylim([-.5 10.5])
axis square

subplot(222)
plot(t,asample(tn,:),'linewidth',2)
hold on
plot(t,gen_fr(tn,:),'.m','linewidth',1,'markersize',12)

plot(ft, fr_hat(tn,:),'-g','linewidth',1.5)
this_spks = t(spikes(tn,:));
if ~isempty(this_spks)
    plot(this_spks,11,'.k','markersize',10)
end
plot(xlim,[0 0],'-k')

xlabel('time from stim on (s)')
legend('a','real rate',  'est rate','spikes','location','northoutside')
title('example trial')


subplot(224)
% imagesc(model_noisy(tn).posterior.pdf',...
%     'x',model_noisy(tn).posterior.T,'y',model_noisy(tn).posterior.avals)


imagesc(model(tn).posterior.pdf','x',model(tn).posterior.T,'y',model(tn).posterior.avals)
hold on
colormap(flipud(bone))
axis xy
plot(t,asample(tn,:),'.-r','linewidth',.1)
ylim([-1 1].*3)
plot(xlim,[0 0],'-k')

%%
model = model_noisy;
switch_t0s = [-.45:.01:.55];
which_switch = 'model';
krndt   = .15;
krn    = my_gauss_krn(krndt,'causal',dt,1);
dt = p.dt;
t0s         = .55:.05:1-lag;
dv_bins   = [-3:.25:3];         % triggers rebinning, but not recompiling
dv_bins = avals_edges   ;%dv_bins   = 20;
alignment   = 'stimstart'
ft = t;
ms = 10;
lw = 2;
compute_neuron = 0;
compute_switch = 0;

for ff = 0:2
    base_fr = 50;
    switch ff
        case 0            
            fr_in = (asample);
            minf = 10;
            k = 10;
            b = 1;
            lc = 'k';
            figure(1); clf
            if compute_switch
                figure(2); clf
            end
            faketuning = @(dv,ft) minf + k*glmval([0 b]',dv(:),'identity');

        case 1
            
            fr_in = sign(asample);
            minf = 0;
            rangef = 10;
            b = .1;            
            lc = 'r';
            faketuning = @(dv,ft) minf + rangef*glmval([0 b]',dv(:),'logit');

        case 2
            fr_in = 50*sign(asample);
            lc = 'm';
            minf = 0;
            rangef = 10;
            b = 1000;  
            faketuning = @(dv,ft) minf + rangef*glmval([0 b]',dv(:),'logit');

    end

    frbins = linspace(min(fr_in(:)),max(fr_in(:)), 50);
    res = dyn_fr_dv_map(cellid, 't0s', t0s, 'model', model, 'lag', lag, ...
        'frates', fr_in, 'ft', ft, ...
        'which_trials',1:length(model),'trialnums',1:length(model),...
        'frbins', frbins, 'alignment', alignment,...
        'krn_width', krn_width, 'krn_type', krn_type,...
        'norm_type', norm_type,'shuffle_trials',0,...
        'n_dv_bins',dv_bins,'end_mask_s',0,...
        'average_fr', 0, 'average_a_vals', 0);
    if compute_switch
        res_switch = dyn_fr_dv_map(excellid, 't0s', switch_t0s, ...
            'frates', fr_in, 'ft', ft,...
            'model', model, 'trialnums',1:length(model),...
            'lag', lag, ...
            'frbins', frbins, 'alignment', alignment,...
            'krn_width', krn_width, 'krn_type', krn_type,...
            'norm_type', norm_type,'shuffle_trials',0,...
            'n_dv_bins',dv_bins,'end_mask_s',0,...
            'which_switch',which_switch);
    end
    if compute_neuron
        
        gen_fr = nan(NT,numel(t));
        for tn = 1:NT;
            lessnoise_sample = asample(tn,:);
            gen_fr(tn,1:length(lessnoise_sample)) = faketuning(lessnoise_sample,tvec);
        end
        spikes = gen_fr * p.dt > rand(size(gen_fr)) ;

        x = 0:krn.bin_size:max_T-krn.bin_size;
        fr_hat = nan(NT, numel(x));
        for tn = 1:NT
            [y x] = spike_filter(0,t(spikes(tn,:)), krn.krn, ...
                'pre',t(1),'post',max_T,'kernel_bin_size', krn.bin_size);
            fr_hat(tn,:) = y;
            ind = find(nanfrates(tn,:),1);
            if ~isempty(ind)
                ind = find(t>cell_ft(ind),1);
                
                fr_hat(tn,ind:end) = nan;
            end
        end
        frbins = linspace(min(fr_hat(:)),max(fr_hat(:)), 50);
        res_neuron = dyn_fr_dv_map(cellid, 't0s', t0s, 'model', model, 'lag', lag, ...
            'frates', fr_hat, 'ft', ft, ...
            'which_trials',1:length(model),'trialnums',1:length(model),...
            'frbins', frbins, 'alignment', alignment,...
            'krn_width', krn_width, 'krn_type', krn_type,...
            'norm_type', norm_type,'shuffle_trials',0,...
            'n_dv_bins',dv_bins,'end_mask_s',0,...
            'average_fr', 0, 'average_a_vals', 0);
        
        
    end
    %%
    
    fh = figure(1);
    
    dvlims = [-5 5]
    
    
    ax = subplot(231);
    fgta_line_plot(res,'linewidth',lw,'ax',ax,...
        'plot_field','fr_given_ta','dvlims',dvlims,'plot_colorbar',0)
    if compute_neuron
        fgta_line_plot(res_neuron,'linewidth',lw,'ax',ax,...
            'plot_field','fr_given_ta','dvlims',dvlims,'plot_colorbar',0)
    end
    title('E[r(a,t)]','interpreter','latex')
    xlabel('time from stim on (s)')
    ylabel('FR')
    
    
    ax_r = subplot(232);
    fgta_line_plot(res,'ax',ax_r,'plot_field','fgta_resid','dvlims',dvlims,...
        'linewidth',lw,'plot_colorbar',0)
    if compute_neuron
        fgta_line_plot(res_neuron,'ax',ax_r,'plot_field','fgta_resid','dvlims',dvlims,...
            'linewidth',lw,'plot_colorbar',0)
    end
    title('E[r(a,t)] - E[r(t)]','interpreter','latex')
    xlabel('time from stim on (s)')
    ylabel('\Delta FR')
    
    ax = subplot(233)
    if ff==0
        plot([-3 3],[-3 3],'color',[1 1 1].*.7)
        hold on
    end
    fgta_plot_tuning(res,'plot_field','fga_resid_tmn', ...
        'errorbar_field','fga_resid_std','ax',ax,'dvlims',dvlims,...
        'linecolor', lc, 'linewidth', 1, 'linestyle', '.-', 'markersize', ms)
    hold on
    if compute_neuron
        fgta_plot_tuning(res_neuron,'plot_field','fga_resid_tmn', ...
            'errorbar_field','fga_resid_std','ax',ax,'dvlims',dvlims,...
            'linecolor', lc, 'linewidth', 1, 'linestyle', '-', 'markersize', ms)
    end
    xlim([-4 4])
    %ylim([-3 3])
    pbaspect(ax, [1 1 1])
    title('E[$\Delta$ r(a)]','interpreter','latex')
    xlabel('accumulated value (a)')
    ylabel('\Delta FR')
    
    
    ax = subplot(236);
    %     if ff==0
    %         plot([-3 3],[-3 3],'--k')
    %         hold on
    %     end
    fgta_plot_tuning(res,'plot_field','rank1_ra_n', ...
        'errorbar_field','','ax',ax,'dvlims',dvlims,...
        'linecolor', lc, 'linewidth',1, 'linestyle', '.-')
    hold on
    if compute_neuron
        fgta_plot_tuning(res_neuron,'plot_field','rank1_ra_n', ...
            'errorbar_field','','ax',ax,'dvlims',dvlims,...
            'linecolor', lc, 'linewidth',1, 'linestyle', '-')
    end
    pbaspect(ax, [1 1 1])
    title('rank 1 $\hat{r}(a)$','interpreter','latex')
    xlabel('accumulated value (a)')
    ylabel('\Delta FR')
    box off
    ax.TickDir = 'out';
    xlim(dvlims)
    %ylim([-1 1].*.5)
    xlim([-4  4])
    
    ax = subplot(235);
    plot(res.t0s, res.rank1_mt_n, '.-', 'linewidth', lw, 'color',lc)
    hold on
    if compute_neuron
        plot(res_neuron.t0s, res_neuron.rank1_mt_n, '-',...
            'linewidth', lw, 'color',lc)
    end
    %pbaspect(ax, [2 1 1])
    title('rank 1 $\hat{m}(t)$','interpreter','latex')
    xlabel('time from stim on (s)')
    ylabel('\Delta FR')
    box off
    ax.TickDir = 'out';
    %     cb = colorbar
    %     drawnow
    %     pos = get(ax,'position')
    %     delete(cb)
    %     set(ax,'position',pos)
    %     ylim([0 5])
    
    ax = subplot(234)
    
    plot(1:length(res.rank_var), res.rank_var, '.-', ...
        'color', lc, 'linewidth', lw)
    hold on
    if compute_neuron
        plot(1:length(res.rank_var), res.rank_var, '-',...
            'color', lc, 'linewidth', lw)
    end
    pbaspect(ax, [1 2 1])
    box off
    ax.TickDir = 'out';
    ylim([.7 1])
    xlabel('rank')
    ylabel('variance explained')
    xlim([1 length(res.rank_var)])
    drawnow
    if compute_switch
        
        fh = figure(2);
        
        dvlims = [-5 5]
        lw = 1;
        
        ax = subplot(231);
        fgta_line_plot(res_switch,'linewidth',lw,'ax',ax,...
            'plot_field','fr_given_ta','dvlims',dvlims,'plot_colorbar',0)
        title('E[r(a,t)]','interpreter','latex')
        xlabel('time from stim on (s)')
        ylabel('FR')
        
        
        ax_r = subplot(232);
        fgta_line_plot(res_switch,'ax',ax_r,'plot_field','fgta_resid','dvlims',dvlims,...
            'plot_colorbar',0)
        title('E[r(a,t)] - E[r(t)]','interpreter','latex')
        xlabel('time from stim on (s)')
        ylabel('\Delta FR')
        
        ax = subplot(233)
        fgta_plot_tuning(res_switch,'plot_field','fga_resid_tmn', ...
            'errorbar_field','fga_resid_std','ax',ax,'dvlims',dvlims,...
            'linecolor',lc)
        hold on
        pbaspect(ax, [1 1 1])
        title('E[$\Delta$ r(a)]','interpreter','latex')
        xlabel('accumulated value (a)')
        ylabel('\Delta FR')
        
        ax = subplot(236);
        fgta_plot_tuning(res_switch,'plot_field','rank1_ra_n', ...
            'errorbar_field','','ax',ax,'dvlims',dvlims,'linecolor',lc)
        hold on
        pbaspect(ax, [1 1 1])
        title('rank 1 $\hat{r}(a)$','interpreter','latex')
        xlabel('accumulated value (a)')
        ylabel('\Delta FR')
        box off
        ax.TickDir = 'out';
        
        ax = subplot(235);
        plot(res_switch.t0s, res_switch.rank1_mt_n, 'linewidth', lw, 'color',lc)
        hold on
        %pbaspect(ax, [2 1 1])
        title('rank 1 $\hat{m}(t)$','interpreter','latex')
        xlabel('time from stim on (s)')
        ylabel('\Delta FR')
        box off
        ax.TickDir = 'out';
        %     cb = colorbar
        %     drawnow
        %     pos = get(ax,'position')
        %     delete(cb)
        %     set(ax,'position',pos)
        
        ax = subplot(234)
        plot(1:length(res_switch.rank_var), res_switch.rank_var, 'color', lc, 'linewidth', lw)
        hold on
        pbaspect(ax, [1 2 1])
        box off
        ax.TickDir = 'out';
        ylim([.7 1])
        xlabel('rank')
        ylabel('variance explained')
        xlim([1 length(res_switch.rank_var)])
    end
end
%%


%%
frbinscell = 0:.25:100;
res_switch = dyn_fr_dv_map(excellid, 't0s', switch_t0s, ...
    'model', model, 'trialnums',1:length(model),...
    'lag', lag, ...
    'frbins', frbinscell, 'alignment', alignment,...
    'krn_width', krn_width, 'krn_type', krn_type,...
    'norm_type', norm_type,'shuffle_trials',0,...
    'n_dv_bins',dv_bins,'end_mask_s',0,...
    'which_switch',which_switch);
lc = 'm'
%%




%%
subplot(212)
plot(t(spikes(tn,:)),5,'.r','markersize',10)
hold on
plot(x,fr_hat(1,:),'m')

%%

A = spikes';
newLen = round(size(A, 1) / 2);
B = reshape(sum(reshape(A, 2, 2, 3)), 2, 3);

%%
krn    = my_gauss_krn(.15,'causal',dwndt,1);

fr_hat  = filter(krn.krn, 1, dwn_spks(1:2,:), [], 2);
%%
offset = ceil(length(krn.krn)/2);
fr_hat = fr_hat(:, 2*offset:end-1);
pre = t(1);
post = t(end);
kernel_bin_size = krn.bin_size;
buffered_pre=pre+offset*kernel_bin_size;
x = -buffered_pre:kernel_bin_size:post;
x = x(offset+1:end-1);

%%

%%









%%
figure(1); clf

imagesc(model(tn).posterior.pdf','x',model(tn).posterior.T,...
    'y', model(tn).posterior.avals)
colormap(flipud(bone))

hold on

plot([0 particle_lessnoise.T],asample,'m')

%%
trial = data(tn);

numsteps = round(trial.T/p.dt);
tvec = p.dt:p.dt:trial.T;
n = 1;
a = zeros(n,numsteps);
[cl, cr]    = make_adapted_cat_clicks(trial.leftbups,...
    trial.rightbups, params(5), params(6));

%slight timing offset on click times because make_click_inputs35 uses qfind, whereas compute_full_trial rounds click times and then bins. The following code is a little messy because it attempts to counteract qfind. I should really just not use make_click_inputs35, but the error is only off by 1 dt, so it doesn't really matter.

[difflr, sumlr] = make_click_inputs35(tvec, trial.leftbups, trial.rightbups,...
    cl, cr);
% Initalize with variance
a(1) = sqrt(params(4)).*randn(1);

% run foward euler
for ii=1:numsteps-1
    a(:,ii+1) = a(:,ii)...
        + p.dt*params(1).*a(:,ii)...
        + difflr(ii) + sqrt(sumlr(ii)*params(3)).*randn(n,1) ...
        + sqrt(params(2)*p.dt).*randn(n,1);
end

plot([0 tvec], a, 'r')
%plot(tvec,difflr,'m')