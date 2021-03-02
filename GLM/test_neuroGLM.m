cd('~/code/neuroGLM')
addpath('~/code/pbups_utils/')
addpath('~/code/phys_utils/')
addpath('/Users/oroville/code/ratter/Analysis/helpers')
addpath(genpath(pwd))
addpath('/Users/oroville/code/ratter/Analysis/Adrian/stat/')
addpath('/Users/oroville/code/ratter/Analysis/Adrian/struct/')
addpath('/Users/oroville/code/ratter/Analysis/Adrian/statistics/')
%addpath(genpath('/Users/oroville/code/ratter/Analysis/Adrian/'))
addpath('/Users/oroville/code/ratter/ExperPort/Analysis')
%%
sd = .15;
type = 'causal';
bin_size = .025;
normalize = 1;
krn = my_gauss_krn(sd,type,bin_size,normalize)
pre = 1.25;
post = 3;
spkwdw = [-pre post];
cellid = 18181;
%cellid = 16857;
%cellid = 17784;
sessid = bdata('select sessid from cells where cellid={S}',cellid);
[ratname, sess_date, pd] = bdata(['select ratname, sessiondate, protocol_data ' ...
    'from sessions where sessid={S}'],sessid);
%%
pd          = pd{1};
peh         = get_peh(sessid);
pdfields    = fieldnames(pd);
ntrials     = length(peh);
for ff = 1:length(pdfields)
    pd.(pdfields{ff}) = pd.(pdfields{ff})(1:ntrials);
end

state_ts    = get_state_ts(peh,pd);
c_ts        = bdata('select ts from spktimes where cellid={S}',cellid);
my_cells.sessid     = sessid;
my_cells.ratname    = ratname;
my_cells.sess_date  = sess_date;
my_cells.mat_file_name = 'neuroGLM_test_temp_db.mat';
my_cells.cellid = cellid;
my_cells.c_ts   = c_ts;
%%
% plot an example psth first

ref = state_ts.stim_start;
pre = 3; post = 2; 
[y, x] = spike_filter(ref, my_cells.c_ts{1}, krn.krn,'kernel_bin_size',krn.bin_size,...
    'pre',pre,'post',post);
[yrast, xrast] = spike_filter(ref, my_cells.c_ts{1}, 1,'kernel_bin_size',.05,'normalize_krn',false,...
    'pre',pre,'post',post);
pokedR = (pd.sides=='r') == pd.hits;
valid_left = ~pd.violations & ~pokedR;
valid_right = ~pd.violations & pokedR;

figure(1234); clf
subplot(221)
imagesc(y,'x',x)
subplot(223)
plot(x,nanmean(y(valid_right,:)),'r')
hold on
plot(x,nanmean(y(valid_left,:)),'b')
plot([0 0],ylim,'k')
drawnow

ref = state_ts.cpoke_end;

cellno = 1;
use_spike_history = 1;

[switch_to_0, switch_to_1, array_data, vec_data] = ...
    get_switches(cellid, 'which_switch','model',...
    'clear_bad_strengths', 1, ...
    'bad_strength', 0, 'fit_line', 1,...
    'exclude_final', exclude_final, 'final_only', final_only);


%%
use_model_switches = 1;
use_generative_switches = 0;
[ex_stats, ex_expt, dm, S, ex_rawData] = fit_glm_to_Cells(my_cells,...
    'useGPU',false,'task','dynamic','cellno',cellno,'nClickBins',1,...
    'use_spike_history',use_spike_history,'pd',pd,'peh',peh,...
    'spikeWindowS',spkwdw, 'use_model_switches',use_model_switches,...
    'use_generative_switches',use_generative_switches,'skip_fit',0);
%%
tn = 201;
array_data(tn).model_switch_to_0;
array_data(tn).switch_to_0;

figure(1); clf
subplot(211)
plot(array_data(tn).left_bups,-1,'.b')
hold on
plot(array_data(tn).right_bups,1,'.r')

ylim([-1 1].*3)
plot(xlim, [0 0], 'k')
plot(array_data(tn).model_T, array_data(tn).model_mean)
if ~isempty(array_data(tn).switch_to_0)
plot(array_data(tn).switch_to_0.*[1 1]' , ylim, 'b')
end
if ~isempty(array_data(tn).switch_to_1)
    plot(array_data(tn).switch_to_1.*[1 1]', ylim, 'r')
end
if ~isempty(array_data(tn).model_switch_to_0)
plot(array_data(tn).model_switch_to_0.*[1 1]', ylim, '--b')
end
if ~isempty(array_data(tn).model_switch_to_1)
    plot(array_data(tn).model_switch_to_1.*[1 1]', ylim, '--r')
end

subplot(212)
plot(ex_rawData.trial(tn).left_clicks,-1,'.b')
hold on
plot(ex_rawData.trial(tn).right_clicks,1,'.r')
ylim([-1 1].*3)
plot(xlim, [0 0], 'k')
%plot(array_data(tn).model_T, array_data(tn).model_mean)
if isfield(ex_rawData.trial(tn),'switch_to_0') & ...
        ~isempty(ex_rawData.trial(tn).switch_to_0)
plot(ex_rawData.trial(tn).switch_to_0.*[1 1]' , ylim, 'b')
end
if isfield(ex_rawData.trial(tn),'switch_to_1') & ...
        ~isempty(ex_rawData.trial(tn).switch_to_1)
    plot(ex_rawData.trial(tn).switch_to_1.*[1 1]', ylim, 'r')
end
if isfield(ex_rawData.trial(tn),'model_to_0') & ...
        ~isempty(ex_rawData.trial(tn).model_to_0)
plot(ex_rawData.trial(tn).model_to_0.*[1 1]', ylim, '--b')
end
if isfield(ex_rawData.trial(tn),'model_to_1') & ...
        ~isempty(ex_rawData.trial(tn).model_to_1)
    plot(ex_rawData.trial(tn).model_to_1.*[1 1]', ylim, '--r')
end

%%
figure(10) 
clf
imagesc(dm.X)
colormap(flipud(bone))


%%
expt = ex_expt;
figure(4568);clf;
drawnow
gammas          = [expt.trial.gamma];

unique_gamma    = unique([expt.trial.gamma]);

[a b bins] = histcounts(abs(gammas),[0 1 2 5]);

[a b bins] = histcounts(abs(gammas),[0 5]);

bins        = bins .* sign(gammas);
unique_bins = unique(bins);
nbins       = length(unique_bins);
cm = colormapRedBlue(nbins/2);
cm(floor(nbins/2)+1,:) = [];
s = [];
hits = [expt.trial.ishit];
base_ind = true(size(hits));
for bb = 1:nbins
    
    tind = find(bins == unique_bins(bb) & base_ind);
    s(1) = subplot(211);
    plotGLM.plotPETH(expt,['sptrain',num2str(cellno)],tind,...
        'color',cm(bb,:),'linestyle','-','linewidth',1);
    hold on
    s(2) = subplot(212);
    plotGLM.plotPETH(expt,S.Yhat,tind,...
        'color',cm(bb,:));
    hold on
    
end
linkaxes(s)
%%
dp.left_color
clf
plotGLM.plotPETH(expt,['sptrain',num2str(cellno)],find(~[expt.trial.pokedR]),...
    'color',dp.left_color,'linestyle',':','linewidth',1);
hold on
plotGLM.plotPETH(expt,['sptrain',num2str(cellno)],find([expt.trial.pokedR]),...
    'color',dp.right_color,'linestyle',':','linewidth',1);

plotGLM.plotPETH(expt,S.Yhat,find(~[expt.trial.pokedR]),...
    'color',dp.left_color,'spike_scale',1000,...
    'plot_errorbar',0,'linewidth',2);

plotGLM.plotPETH(expt,S.Yhat,find([expt.trial.pokedR]),...
    'color',dp.right_color,'spike_scale',1e3,...
    'plot_errorbar',0,'linewidth',2);

%%
fig     = figure(14); clf;

dspec   = ex_stats.dspec;
ws      = ex_stats.ws;
wvars   = ex_stats.wvars;

nCovar = numel(dspec.covar);
for kCov = 1:nCovar
    label = dspec.covar(kCov).label;
    if contains(label,'left') | contains(label,'to_0')
        clr = dp.left_color;
    elseif contains(label,'right') | contains(label,'to_1')
        clr = dp.right_color;
    else
        clr = [.5 .5 .5];
    end
    subplot(4,4, kCov);
    shadedErrorBar(ws.(label).tr/1000, ...
            ws.(label).data, sqrt(wvars.(label).data),{'color',clr});
    
        xlabel(['time from ' strrep(label,'_',' ') ' (s)']);
    hold on
    plot(xlim, [ 0 0],'k')
end
