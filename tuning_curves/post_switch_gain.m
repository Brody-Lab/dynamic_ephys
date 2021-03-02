dp              = set_dyn_path;
cout_auc_file   = fullfile(dp.ephys_summary_dir,'cout_auc.mat');
v               = load(cout_auc_file,'cellids','good_cells');
pop_cellids     = v.cellids(v.good_cells);
pop_cellids     = pop_cellids(1:end-1);

pop_sessids = zeros(size(pop_cellids));

for cc = 1:length(pop_cellids)
    pop_sessids(cc) = bdata(['select sessid from cells where ' ...
        'cellid={S}'],pop_cellids(cc));
end
%%
unique_sessid = unique(pop_sessids);
%%
ops.mask_after_stim_firing    = 1;
ops.mask_stim_firing          = 0;
ops.average_a_vals            = 1;
ops.average_fr                = 1;
ops.min_fr_mod                = 1;
ops.fit_time_sigmoid          = 0;
ops.use_nresidual             = 1;
ops.plot_res                  = 3;

fht = 2.5;
fw  = 1.75*fht;
dp          = set_dyn_path;
direction   = 'backward';
frbins      =  0:.5:100;
krn_type    = 'halfgauss';
end_mask_s  = 0.0;
alignment   = 'stimstart'; %'cpokeout'

lag         = 0.1;            %0.2;
max_dur = 1.975;
switch_t0s = [-.55:.025:.55];
dv_bins = [-5.5:1:5.5];
switch alignment
    case 'stimend'
        t0s         = -1-lag:0.025:-0.1;  %
        t0s         = -max_dur:0.025:-0.1-lag;  %
    case 'stimstart'
        t0s         = 0.1:0.025:max_dur-lag;  %
end
which_switch = 'model';
min_switch_t = 0;
max_switch_t = 3;

fig_type = '-dpng';
norm_type = '';
use_fake = 1;
zscore_frates = 0;
demean_frates = 0;

switch_str = '';
do_save = 0;
do_shuffle_trials = 0;
do_shuffle_fixing_choice = 0;
do_shuffle_switches = 0;

if ~isempty(which_switch)
    t0s_ = switch_t0s;
    switch_str = ['_' which_switch(1:5) 'switch'];
    badtind = t0s_ > -.2 & t0s_ < .20;
    
else
    t0s_ = t0s;
    badtind = t0s_ > Inf;
    
end

dvlims = [];
lw = 2;

if isempty(which_switch)
    xlab = 'time from stim on (s)';
else
    xlab = ['time from ' upper(which_switch) ' switch (s)'];
end
switch mt_field
    case 'rank1_mt_n'
        ylab = {'modulation' '(\Delta FR)'}
    case 'rank1_mt_nm'
        ylab = {'modulation' '(fraction of maximum)'}
        
end
[align_strs, align_args] = dyn_align_LUT;
align_ind = strmatch(alignment,align_strs,'exact');
fht = 2.5;
fw = 1.75*fht;

%%
res_fn = @(cellid, data, model, do_shuffle_switches) dyn_fr_dv_map(cellid, ...
    't0s', t0s_, 'data',data, 'model', model, 'lag', lag, ...
    'alignment', alignment, 'krn_width', krn_width, 'krn_type', krn_type,...
    'norm_type', norm_type, 'shuffle_trials', 0,...
    'n_dv_bins',dv_bins, 'end_mask_s', end_mask_s,...
    'demean_frates', demean_frates, 'zscore_frates', zscore_frates,...
    'which_switch', which_switch, 'shuffle_trials', do_shuffle_trials, ...
    'shuffle_switches', do_shuffle_switches,...
    'min_switch_t',min_switch_t, 'max_switch_t', max_switch_t,...
    'shuffle_trials_fixing_choice', do_shuffle_fixing_choice,...
    'average_fr',1);

%%
profile on

nc = length(pop_cellids);

n_shuffles = 100;

clf;

fprintf('working through %i sessions', length(unique_sessid));
%%
load(['/Users/oroville/projects/pbups_dyn/code/dynamic_ephys/tuning_curves/'...
    'post_switch_gain_workspace.mat'])
%%
counter = 0;

run_shuffles    = 0;
real_diff       = nan(nc,1);

fh = figure(10);

tic
for si = 1:length(unique_sessid)
    toc
    tic
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
        %%
        
        if ii == 1
            model   = get_data_model_p(this_sessid, data.trials.trialnums);
        end
        
        % compute the unshuffled results
        res = res_fn(excellid, data, model, do_shuffle_switches);
        mt_real             = res.(mt_field);
        mt_real(badtind)    = nan;
        real_diff(cc)       = nanmean(mt_real(badn:end)) - ...
            nanmean(mt_real(1:bad0));
        % compute shuffled results
        ra_field = 'rank1_ra_n';
        mt_field = 'rank1_mt_n';
        if run_shuffles
            mt_shuffle = nan(n_shuffles, numel(res.(mt_field)));
            ra_shuffle = nan(n_shuffles, numel(res.(ra_field)));
            tic
            fprintf('shuffle 1...')
            parfor ss = 1:n_shuffles
                if mod(ss,25) == 0
                    fprintf('%i...',ss)
                end
                
                res_shuff = res_fn(excellid, data, model, 1);
                
                ra_shuffle(ss,:) = res_shuff.(ra_field);
                mt_shuffle(ss,:) = res_shuff.(mt_field);
            end
            toc
            if counter == 1
                ptiles = nan(nc,length(res.(mt_field)));
                diff_p = nan(nc,1);
            end
            %%
            shuff_diff    = nanmean(mt_shuffle(:,badn:end),2) - ...
                nanmean(mt_shuffle(:,1:bad0),2);
            ptile           = mean(mt_real' > mt_shuffle);
            ptiles(cc,:)    = ptile;
            diff_p(cc)      = mean(real_diff(cc)>shuff_diff);
            
            clf(fh);
            ax = [];
            %%
            for do_average = 0:1
                %%
                ax(do_average+1) = subplot(1,2,do_average+1);
                this_ax = ax(do_average+1);
                hold(this_ax,'on')
                
                xx = res.t0s;
                xx(badtind) = nan;
                mt_shuff = mt_shuffle;
                mt_shuff(:,badtind) = nan;
                bad0 = find(badtind,1)-1;
                badn = find(badtind,1,'last')+1;
                hold on
                if do_average
                    plot(this_ax, 0, shuff_diff, '.',...
                        'markersize',10,'color', [1 1 1].*.7)
                    plot(this_ax, [-.25 .25], [1 1]*real_diff(cc), '-','color', 'r','linewidth',2)
                    set(gca,'xtick',[])
                else
                    plot(this_ax, xx, mt_shuff, 'color', [1 1 1].*.7, 'linewidth', 1)
                    plot(this_ax, xx([bad0 badn]), mt_shuff(:,[bad0 badn]), ...
                        ':', 'color',[ 1 1 1].*.6,'linewidth', 1)
                    plot(this_ax, xx, mt_real, 'color', [1 1 1].*.0, 'linewidth', lw)
                    yyup    = mt_real;
                    yydwn   = mt_real;
                    yyup(ptile < .95 & ptile > .05)   = nan;
                    yydwn(ptile < .95 & ptile > .05)  = nan;
                    plot(this_ax, xx, yyup, '-', 'color', 'r','linewidth',3)
                    plot(this_ax, xx, yydwn, '-', 'color', 'r','linewidth',3)
                    plot(this_ax, xx([bad0 badn]), mt_real([bad0 badn]), ':', 'color',[ 1 1 1].*0,...
                        'linewidth', lw)
                end
                
                
            end
            
            title('rank 1 $\hat{m}(t)$','interpreter','latex')
            if do_average
                xlabel(this_ax,'')
                ylabel(this_ax,{'increase in mean modulation' 'following switch'})
                set(this_ax,'xtick',[])
            else
                xlabel(this_ax, xlab)
                ylabel(this_ax, ylab)
            end
            box off
            set(this_ax, 'YGrid', 'on');
            %             linkaxes(ax,'x')
            %             xlim(ax(1),res.t0s([1 end])+[-.1 .1])
            %             ylim(this_ax,[min(min(ylim,0)) max(ylim)])
        end
    end
end
%%
figure(2); clf
sig = diff_p>.95;
plot(real_diff,diff_p,'ob','markerfacecolor',[.85 .88 1])
hold on
scatter(real_diff(sig),diff_p(sig),'ob','markerfacecolor',[1 .25 .25])
sig = diff_p<.05;
scatter(real_diff(sig),diff_p(sig),'ob','markerfacecolor',[.25 1 .25])

plot([0 0],ylim,'k--')
ylabel('%tile')
xlabel('mean modulation difference (post-pre)')
box off
grid on
axis tight
axis square

%%
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
