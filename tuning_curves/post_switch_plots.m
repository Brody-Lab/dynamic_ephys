
sp = struct('which_switch', 'model',...
            'clear_bad_strengths', 1, 'bad_strength', 0, 'fit_line', 1,...
            't_buffers',[.2 .2], 'min_pre_dur', 0, 'min_post_dur', 0,...
            'min_switch_t',0,'max_switch_t',Inf,...
            'exclude_final',0,'final_only',0,'model_smooth_wdw',100);
this_stadir = get_sta_dirname(sp);
load(fullfile(this_stadir,   'post_switch_gain_workspace.mat'))
%%
do_print = 0;
fht = 2.5;
fw = 1.75*fht;
alpha = .05;
fh = figure(2); clf
set(fh,'position',[10 10 fw fht])
this_ax = axes();
edges = [0:alpha/2:1]

% histogram(this_ax,diff_p,edges,'facecolor',[1 1 1])
% hold on
% histogram(this_ax,diff_p(diff_p>.95),edges,'facecolor',[1 1 1].*0)
ncells = sum(~isnan(diff_p));
h = histogram(this_ax,diff_p,edges,'facecolor',[1 1 1])
hold on
histogram(this_ax,diff_p(diff_p>(1-alpha/2)),h.BinEdges,'facecolor',[1 1 1].*0)
histogram(this_ax,diff_p(diff_p<(alpha/2)),h.BinEdges,'facecolor',[1 1 1].*0)
plot([0 1],[1 1].*	,'r--')
ylabel('# cells')
xlabel(this_ax,{'modulation increase percentile rel. shuffle' })
title('')
hl=legend('all cells','significant')
set(hl,'location','eastoutside')
box off
axis tight
title('switch histogram','fontweight','normal')
%plot([1 1].*mean(real_diff),ylim,'k')
if do_print
    print(fh,fullfile(dp.fig_dir,...
        ['gain_ptile_hist']),'-dsvg','-painters');
end

    
%%
hi_sig = diff_p>(1-alpha/2);
lo_sig = diff_p<alpha/2;

fh = figure(3); clf
set(fh,'position',[10 5 fw fht],'paperposition',[0 0 fw fht],...
    'papersize',[fw fht])
this_ax = axes();
edges = [0:.1:1]

% histogram(this_ax,diff_p,edges,'facecolor',[1 1 1])
% hold on
% histogram(this_ax,diff_p(diff_p>.95),edges,'facecolor',[1 1 1].*0)
edges = [-10:.75:10]
h = histogram(this_ax,real_diff,edges,'displaystyle','stairs','edgecolor','k')
h = histogram(this_ax,real_diff,edges,'facecolor',[1 1 1].*.95,...
                                        'edgecolor',[1 1 1].*.45)
hold on
histogram(this_ax,real_diff(hi_sig | lo_sig),h.BinEdges,...
      'edgecolor',[1 1 1].*.45,...
      'facecolor',[1 1 1].*.25)
% histogram(this_ax,real_diff(diff_p>.95 | diff_p<.05),h.BinEdges,...
%      'displaystyle','stairs', 'edgecolor','k')

%title('switch histogram','fontweight','normal')

ylabel('# cells')
xlabel(this_ax,{'modulation increase (post-pre)' })
hl=legend('all cells','significant')
set(hl,'location','eastoutside')
box off

axis tight
box(hl, 'off')
set(this_ax,'position',get(ex_ax,'position'))
pbaspect(this_ax,[1 1 1])

%plot([1 1].*mean(real_diff),ylim,'k')
print(fh,fullfile(dp.fig_dir,...
        ['gain_mag_hist']),'-dsvg','-painters');
    
%%


%%
excellid = 18181;
data = dyn_cell_packager(excellid);
model   = get_data_model_p(data, data.trials.trialnums);

res = res_fn(excellid, data, model, 0);
mt_shuffle = nan(n_shuffles, numel(res.(mt_field)));
ra_shuffle = nan(n_shuffles, numel(res.(ra_field)));

for ss = 1:n_shuffles
    if mod(ss,20) == 0 
        fprintf('shuffle %i',ss)
    end
    res_shuff = res_fn(excellid, data, model, 1);
    
    ra_shuffle(ss,:) = res_shuff.(ra_field);
    mt_shuffle(ss,:) = res_shuff.(mt_field);
end
%%

mt_real = res.(mt_field);
ex_shuff_diff    = nanmean(mt_shuffle(:,badn:end),2) - ...
    nanmean(mt_shuffle(:,1:bad0),2);
ex_real_diff  = nanmean(mt_real(badn:end)) - ...
    nanmean(mt_real(1:bad0));

%%
fh = figure(1); clf
set(fh,'position',[10 10 fw fht])
this_ax = axes();
c = .25
plot(this_ax, c*rand(size(ex_shuff_diff)) - c/2,ex_shuff_diff, '.',...
    'markersize',10,'color', [1 1 1].*.7)
hold on
plot(this_ax, [-.25 .25], [1 1]*ex_real_diff, '-','color', 'k','linewidth',2)
set(this_ax,'xtick',[])
ylabel(this_ax,{'modulation increase' '(post-pre)'})
box off
pbaspect(this_ax,[.2 1 1])
ylim([-3 7])
text(c*1.5, 5, 'data','color','k','fontsize',dp.fsz)
text(c*1.5, 4, 'shuffles','color',[1 1 1].*.6,'fontsize',dp.fsz)
print(fh,fullfile(dp.fig_dir,...
        ['example_gain_shuffle_points']),'-dsvg','-painters');
%%
fh = figure(1); clf
set(fh,'position',[10 5 fw fht],'paperposition',[0 0 fw fht],...
    'papersize',[fw fht])
clf 
ex_ax = axes();
plot([1 1]*ex_real_diff,ylim,'k')
hold on
histogram(ex_ax,ex_shuff_diff,10,'facecolor',[1 1 1].*.7,...
    'edgecolor',[1 1 1].*.7)
xlim([-2.15 7])
pbaspect(ex_ax,[1 1 1])
plot([1 1]*ex_real_diff,ylim,'k')
xlabel(ex_ax,{'modulation increase (post-pre)'})
ylabel('# tests')
box(ex_ax,'off')
hl=legend('data','shuffle    ')
set(hl,'location','eastoutside')
print(fh,fullfile(dp.fig_dir,...
        ['example_gain_shuffle_hist']),'-dsvg','-painters');
%%
%%


%%
switch_ratnames = {};
%%

for cc = 1:length(pop_sessids)
    this_rat = bdata('select ratname from sessions where sessid={S}',pop_sessids(cc));
    switch_ratnames{cc} =   this_rat{1};
end
%%
rats=unique(switch_ratnames,'stable')
b=cellfun(@(x) sum(ismember(switch_ratnames,x)),rats,'un',0)


rats=unique(switch_ratnames(diff_p>.95),'stable')
%%
cell_list = dyn_cells_db;
select_str = 'normmean > 0'
ratnames = extracting(cell_list, 'ratname', select_str);
prefp = cell2mat(extracting(cell_list, 'prefp', select_str));
normmean = cell2mat(extracting(cell_list, 'normmean', select_str));
sig = prefp < .05;
nm1_cells = normmean > 1;

rats = unique(ratnames,'stable');
%% Rat ephys table
rats = dp.ratlist';
ratnums = cellfun(@(x) str2num(x(2:end)),rats);   
ncells = arrayfun(@(x) sum(ismember(ratnames,x)),rats);
nsig = arrayfun(@(x) sum(ismember(ratnames(sig),x)),rats);
nm1 = arrayfun(@(x) sum(ismember(ratnames(nm1_cells),x)),rats);
nm1xsig = arrayfun(@(x) sum(ismember(ratnames(sig&nm1_cells),x)),rats);

%siggain_rat = arrayfun(@(x) sum(ismember(switch_ratnames(diff_p>.95),x)),rats)';

T = table(rats,ncells,nm1,nsig,nm1xsig);
sortrows(T,'ncells','descend')
%%

for si = 1:length(unique_sessid)
    fprintf('session %i...',si)
    this_sessid = unique_sessid(si);
    ind     = find(pop_sessids == unique_sessid(si));
    for ii = 1:length(ind)
        cc          = ind(ii);
        counter     = counter + 1;
        excellid    = pop_cellids(cc);
        fprintf('\ncell %i...',counter)
        %%
        data = dyn_cell_packager(excellid);
        switch_ratnames{cc} =   data.ratname;
        
        this_rat = bdata('select ratname from sessions where sessid={S}',pop_sessids(cc));
        
        switch_ratnames{cc} =   this_rat{1};
    end
end






%% FRM vs post gain boost
t0s = .01:.025:1.25;
t0_res_fn = @(cellid, data, model) dyn_fr_dv_map(cellid, ...
    'data', data, 'model', model, ...
    't0s', t0s,  'lag', lag, ...
    'alignment', alignment, 'krn_width', krn_width, 'krn_type', krn_type,...
    'norm_type', norm_type, 'shuffle_trials', 0,...
    'n_dv_bins',dv_bins, 'end_mask_s', end_mask_s,...
    'demean_frates', demean_frates, 'zscore_frates', zscore_frates,...
    'which_switch', which_switch, 'shuffle_trials', do_shuffle_trials, ...
    'shuffle_switches', do_shuffle_switches,...
    'min_switch_t',min_switch_t, 'max_switch_t', max_switch_t,...
    'shuffle_trials_fixing_choice', do_shuffle_fixing_choice,...
    'average_fr',1);
counter = 0;
frm = nan(nc, 1);
for si = 1:length(unique_sessid)
    fprintf('session %i...',si)
    this_sessid = unique_sessid(si);
    ind     = find(pop_sessids == unique_sessid(si));
    for ii = 1:length(ind)
        counter = counter + 1;
        cc          = ind(ii);
        excellid    = pop_cellids(cc);
        fprintf('\ncell %i...',counter)
        %%
        data = dyn_cell_packager(excellid);
        %%
        if ii == 1
            model   = get_data_model_p(this_sessid, data.trials.trialnums);
        end
        
        % compute the unshuffled results
        res = t0_res_fn(excellid, data, model);
        frm(cc) = res.fr_mod;
    end
end
%%
figure(1); clf
sig = diff_p>.95;
b = glmfit(frm(:),diff_p(:),'normal');
xx = [0:.1: 15];
% plot(xx, glmval(b,xx,'identity'),'--','color',[1 1 1].*.5)
% b = glmfit(frm(:),diff_p(sig),'normal');
% xx = [0:.1: 15];
plot(xx, glmval(b,xx,'identity'),'--','color',[1 .5 .5].*.5)
hold on
plot(frm,diff_p,'ob','markerfacecolor',[.85 .88 1])

scatter(frm(sig),diff_p(sig),'ob','markerfacecolor',[1 .25 .25])
box off
ylabel('%tile')
xlabel('mean modulation (spks/s)')
