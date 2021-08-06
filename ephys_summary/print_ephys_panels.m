close all
clear all
dp = set_dyn_path(1)
%%
% produce a stim on aligned plot and a center out aligned plot
align_strs      = dyn_align_LUT;
cout_align_ind  = find(ismember(align_strs,'cpokeout'));
cin_align_ind   = find(ismember(align_strs,'stimstart-cout-mask'));
assert(all([~isempty(cout_align_ind),~isempty(cin_align_ind)]));

% turn long_trials_only on if you want to only use the long trials for the
% figures
long_trial_dur = 1;
long_trials_only = false;

%% get full cell_list; *all* units recorded in project
cell_list = dyn_cells_db;   % That has to be run once to create cell_list
normmnth = 1;
prefpth  = .05;

%%
fht = 2.5;
fw = 1.75*fht;
fsz = 10.5;
set(0, 'defaultaxesfontsize',fsz);
set(0,'defaultaxeslinewidth',1)
xlim_on = [-.55 1.5];
xlim_off = [-1.55 .5];
xlim_on = [-.55 1.25];
xlim_off = [-1.25 .55];
%%
[~, align_args]      = dyn_align_LUT;
cellid = 18181;
[ad, vd] = package_dyn_phys(cellid);
i = 2;
[sp_counts{i}, sp_count_T{i}] = calc_sp_counts(cellid, ...
        'array_data', ad, 'vec_data', vd, align_args{i}{:});

%% load up spikes and make psths for each cell
recompute = 0;
rind = 3;
lind = 2;
    
all_psth_fn = fullfile(dp.spikes_dir, 'all_cells_psth.mat');
if exist(all_psth_fn, 'file') & ~recompute
    fprintf('loading file...\n')
    load(all_psth_fn);
else
    select_str = 'normmean > 0' ;
    ratnames = cell2mat(extracting(cell_list, 'ratname', select_str));
    cellids = cell2mat(extracting(cell_list, 'cellid', select_str));
    sessids = cell2mat(extracting(cell_list, 'sessid', select_str));

    psr = cell2mat(extracting(cell_list, 'prefsideright', select_str));
    prefp = cell2mat(extracting(cell_list, 'prefp', select_str));
    normmean = cell2mat(extracting(cell_list, 'normmean', select_str));
    ncells = size(cellids,1);
    disp(['loaded ' num2str(ncells) ' cells']);
    
    for cc = 1:ncells
        
        if mod(cc,25)==0
            fprintf([num2str(cc) '...'])
        end
        d=dyn_cell_packager(cellids(cc));
        cin_frates = d.frate{cin_align_ind};
        cout_frates = d.frate{cout_align_ind};
        
        % now that we know how many timepoints we have, initialize variables
        if cc == 1
            cin_t = d.frate_t{cin_align_ind};
            cout_t = d.frate_t{cout_align_ind};
            cin_binsz = diff(cin_t([1 2]));
            cout_binsz = diff(cout_t([1 2]));
            cin_ntp = length(cin_t);
            cout_ntp = length(cout_t);
            
            % 3rd dimension is both, l, r, both
            % 4th dimension is all data, split A, split B
            cin_psth = nan(length(cellids),cin_ntp,3,3);
            cout_psth = nan(length(cellids),cout_ntp,3,3);
            cin_psth_hit = nan(length(cellids),cin_ntp,3,3); % 3rd dimension is both, l, r, both
            cout_psth_hit = nan(length(cellids),cout_ntp,3,3);
            cin_psth_err = nan(length(cellids),cin_ntp,3,3); % 3rd dimension is both, l, r, both
            cout_psth_err = nan(length(cellids),cout_ntp,3,3);
            
        end
        
        poke_r = d.trials.rat_dir==1;
        poke_l = d.trials.rat_dir==-1;
        hit = d.trials.hit==1;
        err = hit == 0;
        stim_dur = d.trials.cpoke_end - d.trials.stim_start;
        good = true(size(stim_dur));
        if long_trials_only
            good = stim_dur > long_trial_dur;
        else
            good_ind = find(good);
        end
        split_ind = randperm(length(good_ind));
        splitAind = split_ind(1:floor(end/2));
        splitBind = split_ind(ceil(end/2):end);
        
        for xx = 1:3
            if xx == 2
                good = false(size(stim_dur));
                good(splitAind) = true;
            elseif xx==3
                good = false(size(stim_dur));
                good(splitBind) = true;
            end
            % hits & errors combined
            cin_psth(cc,:,1,xx)    = nanmean(cin_frates(good,:));
            cin_psth(cc,:,lind,xx)    = nanmean(cin_frates(good&poke_l,:));
            cin_psth(cc,:,rind,xx)    = nanmean(cin_frates(good&poke_r,:));
            
            cout_psth(cc,:,1,xx)   = nanmean(cout_frates(good,:));
            cout_psth(cc,:,lind,xx)   = nanmean(cout_frates(good&poke_l,:));
            cout_psth(cc,:,rind,xx)   = nanmean(cout_frates(good&poke_r,:));
            
            % hits only
            cin_psth_hit(cc,:,1,xx)    = nanmean(cin_frates(hit&good,:));
            cin_psth_hit(cc,:,lind,xx)    = nanmean(cin_frates(hit&good&poke_l,:));
            cin_psth_hit(cc,:,rind,xx)    = nanmean(cin_frates(hit&good&poke_r,:));
            
            cout_psth_hit(cc,:,1,xx)   = nanmean(cout_frates(hit&good,:));
            cout_psth_hit(cc,:,lind,xx)   = nanmean(cout_frates(hit&good&poke_l,:));
            cout_psth_hit(cc,:,rind,xx)   = nanmean(cout_frates(hit&good&poke_r,:));
            
            % errors only
            cin_psth_err(cc,:,1,xx)    = nanmean(cin_frates(err&good,:));
            cin_psth_err(cc,:,lind,xx)    = nanmean(cin_frates(err&good&poke_l,:));
            cin_psth_err(cc,:,rind,xx)    = nanmean(cin_frates(err&good&poke_r,:));
            
            cout_psth_err(cc,:,1,xx)   = nanmean(cout_frates(err&good,:));
            cout_psth_err(cc,:,lind,xx)   = nanmean(cout_frates(err&good&poke_l,:));
            cout_psth_err(cc,:,rind,xx)   = nanmean(cout_frates(err&good&poke_r,:));
        end
        
    end
    save(all_psth_fn, 'cin_t', 'cout_t',...
        'cin_psth', 'cout_psth',...
        'cin_psth_hit', 'cout_psth_hit',...
        'cin_psth_err', 'cout_psth_err',...
        'ratnames', 'cellids', 'sessids', 'psr', 'prefp', 'normmean', 'ncells');
end
%%
fprintf('number of cells: %i \n number of sessions: %i',...
    length(cellids), length(unique(sessids)))
%%
rats = dp.ratlist';

ncells = nan(size(rats));
nsess = nan(size(rats));
for rr = 1:length(rats)
    idx = sum(ratnames == rats{rr},2)==4;
    ncells(rr) = sum(idx);
    nsess(rr) = length(unique(sessids(idx)));
end
fprintf('\nmin/max cells: %i/%i \nmin/max sessions: %i/%i\n',...
    min(ncells),max(ncells),min(nsess),max(nsess))
fprintf('\nmean cells: %.2f \nmean sessions: %.2f\n',...
    mean(ncells),mean(nsess))
table(rats,ncells,nsess)
%% plot example cells in panel B
ppos = [8 10 fw fht ]
cellid = 18181;
[fh, ax] = example_cell_psth('cells',cellid,...
    'cintrange',xlim_on,'couttrange',xlim_off,'coutstr','stimend-nomask','fig_num',2)

ylim(ax,[14 60])
xlim(ax(1),xlim_on)
xlim(ax(2),xlim_off)
ax(1).YTick = [15:15:60];
ax(2).YTick = [15:15:60];
ax(2).YTickLabel = {};

set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],...
    'papersize',ppos([3 4]))
%pbaspect(ax(1),[1 .8 1])
%pbaspect(ax(2),[1 .8 1])

print(fh, fullfile(dp.psth_fig_dir, ['cell_' num2str(cellid) ]),...
    '-dsvg', '-painters')

%%
cellid = 16857;
[fh ax] = example_cell_psth('cells',cellid)%,'coutstr','stimend-no')
% [fh ax] = example_cell_psth('cells',cellid,...
%      'cintrange',xlim_on,'couttrange',xlim_off,'coutstr','stimend-nomask','fig_num',2)

ylim(ax,[0 20])
xlim(ax(1),xlim_on)
xlim(ax(2),xlim_off)
% ax(1).YTick = [15:15:60];
% ax(2).YTick = [15:15:60];
ax(2).YTickLabel = {};
set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]))
pbaspect(ax(1),[1 .8 1])
pbaspect(ax(2),[1 .8 1])

print(fh, fullfile(dp.psth_fig_dir, ['cell_' num2str(cellid) ]),...
    '-dsvg', '-painters')

%%
cellid = 17784;
[fh ax] = example_cell_psth('cells',cellid)
% [fh, ax] = example_cell_psth('cells',cellid,...
%     'cintrange',xlim_on,'couttrange',xlim_off,'coutstr','stimend-nomask','fig_num',2)

ylim(ax,[0 40])
xlim(ax(1),xlim_on)
xlim(ax(2),xlim_off)
% ax(1).YTick = [15:15:60];
% ax(2).YTick = [15:15:60];
ax(2).YTickLabel = {};
pbaspect(ax(1),[1 .8 1])
pbaspect(ax(2),[1 .8 1])
set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]))
print(fh, fullfile(dp.psth_fig_dir, ['cell_' num2str(cellid) ]),...
    '-dsvg', '-painters')




%% plot center in and out aligned plots sorted
% for now, I'm abandoning the sequences with these data, because the rank
% correlations between sequential orderings in splits of the data are quite
% different and I don't know what to make of that

% determine whether to plot the sequence using the sorting data or using
% the heldout data
plot_sorting_data = true;
% decide which timepoints to use for the plot
cin_ind = cin_t>-.25 & cin_t<2;%5;
cout_ind = cout_t>-.25;
save_2name = 'sequence_plot';

good_cells = mean(cin_psth_hit(:,cin_ind,1,1),2) > 2;
pref_r = nanmean(cin_psth_hit(:,:,rind,1),2) > nanmean(cin_psth_hit(:,:,lind,1),2);

cin_psth_hit_pref           = vertcat(cin_psth_hit(good_cells&pref_r,:,rind,1), ...
    cin_psth_hit(good_cells&~pref_r,:,lind,1));
cout_psth_hit_pref          = vertcat(cout_psth_hit(good_cells&pref_r,:,rind,1),...
    cout_psth_hit(good_cells&~pref_r,:,lind,1));
cin_psth_hit_pref_xval1     = vertcat(cin_psth_hit(good_cells&pref_r,:,2,2), ...
    cin_psth_hit(good_cells&~pref_r,:,lind,2));
cout_psth_hit_pref_xval1    = vertcat(cout_psth_hit(good_cells&pref_r,:,2,2), ...
    cout_psth_hit(good_cells&~pref_r,:,lind,2));
cin_psth_hit_pref_xval2     = vertcat(cin_psth_hit(good_cells&pref_r,:,2,3),...
    cin_psth_hit(good_cells&~pref_r,:,lind,3));
cout_psth_hit_pref_xval2    = vertcat(cout_psth_hit(good_cells&pref_r,:,2,3),...
    cout_psth_hit(good_cells&~pref_r,:,lind,3));

combo_psthA = [cin_psth_hit_pref_xval1(:,cin_ind) cout_psth_hit_pref_xval1(:,cout_ind)];
combo_psthB = [cin_psth_hit_pref_xval2(:,cin_ind) cout_psth_hit_pref_xval2(:,cout_ind)];
%combo_psthA = [cin_psth_hit_pref_xval1(:,:) cout_psth_hit_pref_xval1(:,:)];
%combo_psthB = [cin_psth_hit_pref_xval2(:,:) cout_psth_hit_pref_xval2(:,:)];

[combo_sortedA, sortA, sortAt] = sort_by_peak(combo_psthA);
[combo_sortedB, sortB] = sort_by_peak(combo_psthB);
if plot_sorting_data
    psth_all = [cin_psth_hit_pref_xval1(:,cin_ind) cout_psth_hit_pref_xval1(:,cout_ind)];
    %    psth_all = [cin_psth_hit_pref_xval1(:,:) cout_psth_hit_pref_xval1(:,:)];
    for i=1:size(psth_all,1)
        psth_all(i,:) = norm_by_peak(psth_all(i,:));
    end
    
    psth_cin = psth_all(:,1:sum(cin_ind));
    psth_cout = psth_all(:,sum(cin_ind)+1:end);
    %    psthL = psth_all(:,1:120);
    %    psthR = psth_all(:,121:end);
    %    psthL = psthL(:,cin_ind);
    %    psthR = psthR(:,cout_ind);
    
else
    psth_cin = cin_psth_hit_pref_xval2(:,cin_ind);
    psth_cout = cout_psth_hit_pref_xval2(:,cout_ind);
end

tA = cin_t(cin_ind);
tB = cout_t(cout_ind);


%% plot chronometric population psth 
mn_fr       = nanmean(cin_psth(:,:,1,1),2);
mn_fr       = nanmean(cin_psth(:, cin_t > 0, 1, 1),2);
active_cells = mn_fr > normmnth;
sig_cells = prefp < prefpth;
good_cells  = active_cells & sig_cells;
fprintf(['n good cells = ' num2str(sum(good_cells))])
fprintf('\n%.1f %% (%i/%i) of active cells are signficant',...
    100*sum(good_cells)/sum(active_cells),sum(sig_cells),sum(active_cells))

edges = [-2 -.5 -.25  0  .25 .5  2];
pref_color  = [.8 .25 .8];
npref_color = [.8 .65 .25];
[fh, ax] = example_cell_psth('separate_hits', 0, 'min_t', 0, ...
    'cells', cellids(good_cells), 'meta', 1, 'type', 'chrono', ...
    'edges', edges, 'norm', 'onset', 'top_color', pref_color, ...
    'bot_color', npref_color,'fig_num', 1);

xlim(ax(1),xlim_on)
xlim(ax(2),xlim_off)
txt = text(ax(1), .5, 2, ['n = ' num2str(sum(good_cells))])
ax(2).YTickLabel = {};
ylim(ax(1),[.9 2.5])
pbaspect(ax(1),[1 1 1])
pbaspect(ax(2),[1 1 1])

set(fh,'position',ppos,...
    'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]))
groupname = sprintf('cells_normmnth_%i_prefpth_%.2f.svg',normmnth,prefpth);
%%
print(fh, fullfile(dp.psth_fig_dir, groupname),...
    '-dsvg', '-painters')
%% load auc for all neurons
cout_auc_file = fullfile(dp.ephys_summary_dir,'cout_auc.mat');
recompute = 0;
nboot = 250;

if ~exist(cout_auc_file) | recompute
    dp_cin = nan(size(cin_psth(:,:,1,1)));
    dp_cout = nan(size(cout_psth(:,:,1,1)));
    
    good_cell_ix = 1:ncells;
    ngood = length(good_cell_ix);
    
    n_coutt = length(cout_t);
    
    cout_auc = nan(ngood,n_coutt);
    cout_p   = nan(ngood,n_coutt);
    cout_ci  = nan(ngood,n_coutt,2);
    fprintf('working on cell...')
    tic
    for cc = 1:ncells
        if mod(cc,5)==0
            toc;
            fprintf('%i of %i...',cc,ncells);
            tic;
        end
        this_id = cellids(good_cell_ix(cc));
        
        d           = dyn_cell_packager(this_id);
        cin_frates  = d.frate{cin_align_ind};
        cout_frates = d.frate{cout_align_ind};
        go_r        = d.trials.rat_dir==1;
        go_l        = d.trials.rat_dir==-1;
        psth_r      = cout_frates(go_r,:);
        psth_l      = cout_frates(go_l,:);
        
        parfor tt = 1:n_coutt
            [cout_auc(cc,tt), cout_p(cc,tt), cout_ci(cc,tt,:)] = ...
                bootroc(psth_r(:,tt),psth_l(:,tt),nboot);
        end
    end
    save(cout_auc_file,'sessids','cellids','good_cells','cout_t','cout_auc','cout_p','cout_ci','nboot');

else
    load(cout_auc_file);
end

%%
%for 

%% plot sorted auc
%box(s(ii),'off')
    %pbaspect(s(ii),[1 .8 1])
n_consec_sig = 8;
alpha = .05;
x0 = xlim_off(1);
gnames = {'all' 'good' };
fh = figure(10); clf;
fht_ = 1.3*fht;
set(fh,'position',[10 10 fw fht_],'paperposition',[0 0 fw fht_],...
    'papersize',[fw fht_])
for gg = 1:length(gnames)
    %
    subplot(1,2,gg)

    gn = gnames{gg};
    switch gn
        case 'all'
            groupname = 'mnfr_1';
            ind = mn_fr > 1;
        case 'good'
            ind = good_cells;
            groupname = sprintf('cells_normmnth_%i_prefpth_%.2f.svg',normmnth,prefpth);
    end
    
    sigR  = cout_p(ind,cout_t > x0) > 1-alpha/2;
    sigL  = cout_p(ind,cout_t > x0) < alpha/2;
    
    sigsumL = movsum(sigL, n_consec_sig, 2);
    sigsumR = movsum(sigR, n_consec_sig, 2);
    
    sigsumnp  = (sigsumL == n_consec_sig) | (sigsumR == n_consec_sig);
    sigsumnp(:,end) = ones(size(sigsumnp,1),1).*n_consec_sig;
    these_cells = find(ind);
    [sorted_heatplot, sort_id] = sort_by_peak(cout_auc(ind,:),sigsumnp);
    these_cells_sorted = these_cells(sort_id);
    imagesc(sorted_heatplot,'x',cout_t,'Alphadata',~isnan(sorted_heatplot))
    set(gca,'color',[1 1 1].*.85)
    colormap(colormapRedBlue.^.6)
    axpos = get(gca,'position')
    switch gg
        case 1
            title({'all active cells'},'fontweight','normal')
        case 2
            title({'side-selective cells'},'fontweight','normal')
    end
    drawnow
    ylabel({'cell # (sorted by latency)'})

    if gg == 2
        cb = colorbar
        title(cb, 'AUC')
        ylabel('')
        cb.Position = cb.Position + [0.125 0.125 .02 -.5];
    end
    caxis([0.25 .75])
    xlim(xlim_off)
    hold on
    line([0 0],ylim,'color','k')
    box off
    set(gca, 'position', axpos)
    cm = flipud(colormapLinear(dp.left_color,50));
    cm = [cm; colormapLinear(dp.right_color,50)];
    colormap(cm)
    xlabel('time from movement (s)')
   
    print(fh,fullfile(dp.fig_dir,['coutauc' groupname]),'-dsvg','-painters')
end


%% NOT USED - plot PREF/NONPREF pop average PSTH for panel D
mn_fr       = nanmean(cin_psth_hit(:,:,1,1),2);

good_cint   = cin_t >= xlim_on(1) & cin_t <= xlim_on(2);
good_coutt  = cout_t >= xlim_off(1) & cout_t <= xlim_off(2);

switch 0
    case 0
        normnth = 1;
        prefpth = .05;
        good_cells  = mn_fr > 1 & prefp < .05;
        pref_r = nanmean(cin_psth_hit(:,:,rind,1),2) > nanmean(cin_psth_hit(:,:,lind,1),2);
        pref_l = ~pref_r;
    case 1
        normnth = 5;
        prefpth = .01;
        pref_r = psr == 1;
        pref_l = 0   == psr;
    case 2
        normnth = 1;
        prefpth = .05;
        pref_r = psr == 1;
        pref_l = 0   == psr;
end

good_cells  = mn_fr > normnth & prefp < prefpth;
pref_r = psr == 1;
pref_l = 0   == psr;

fprintf(['n good cells = ' num2str(sum(good_cells))])

cin_pref_psth = [cin_psth_hit(pref_r,good_cint,rind,1); ...
    cin_psth_hit(pref_l,good_cint,lind,1)];
cin_npref_psth = [cin_psth_hit(pref_r,good_cint,lind,1); ...
    cin_psth_hit(pref_l,good_cint,rind,1)];

cout_pref_psth = [cout_psth_hit(pref_r,good_coutt,rind,1); ...
    cout_psth_hit(pref_l,good_coutt,lind,1)];
cout_npref_psth = [cout_psth_hit(pref_r,good_coutt,lind,1); ...
    cout_psth_hit(pref_l,good_coutt,rind,1)];

cin_pref_psth_err = [cin_psth_err(pref_r,good_cint,rind,1); ...
    cin_psth_err(pref_l,good_cint,lind,1)];
cin_npref_psth_err = [cin_psth_err(pref_r,good_cint,lind,1); ...
    cin_psth_err(pref_l,good_cint,rind,1)];

cout_pref_psth_err = [cout_psth_err(pref_r,good_coutt,rind,1); ...
    cout_psth_err(pref_l,good_coutt,lind,1)];
cout_npref_psth_err = [cout_psth_err(pref_r,good_coutt,lind,1); ...
    cout_psth_err(pref_l,good_coutt,rind,1)];

switch 0
    case 0
        mn_fr_t0 = 1;
    case 1
        mn_fr_t0 = max([cin_pref_psth cout_pref_psth],[],2);
    case 2
        mn_fr_t0 = normmean([find(pref_r); find(pref_l)]);
end
%good_cells = [cellids(pref_r); cellids(pref_l)] == 18181;


cin_pref_psth = cin_pref_psth./mn_fr_t0;
cin_npref_psth = cin_npref_psth./mn_fr_t0;
cout_pref_psth = cout_pref_psth./mn_fr_t0;
cout_npref_psth = cout_npref_psth./mn_fr_t0;

cin_pref_psth_err = cin_pref_psth_err./mn_fr_t0;
cin_npref_psth_err = cin_npref_psth_err./mn_fr_t0;
cout_pref_psth_err = cout_pref_psth_err./mn_fr_t0;
cout_npref_psth_err = cout_npref_psth_err./mn_fr_t0;

cin_pref_psth_good = cin_pref_psth(good_cells,:);
cin_pref_psth_mn = nanmean(cin_pref_psth_good,1);
cin_npref_psth_good = cin_npref_psth(good_cells,:);
cin_npref_psth_mn = nanmean(cin_npref_psth_good,1);

cout_pref_psth_good = cout_pref_psth(good_cells,:);
cout_pref_psth_mn = nanmean(cout_pref_psth_good,1);
cout_npref_psth_good = cout_npref_psth(good_cells,:);
cout_npref_psth_mn = nanmean(cout_npref_psth_good,1);

cin_pref_psth_good_err = nanmean(cin_pref_psth_err(good_cells,:),1);
cin_npref_psth_good_err = nanmean(cin_npref_psth_err(good_cells,:),1);

cout_pref_psth_good_err = nanmean(cout_pref_psth_err(good_cells,:),1);
cout_npref_psth_good_err = nanmean(cout_npref_psth_err(good_cells,:),1);

psths = [cin_pref_psth_mn cout_pref_psth_mn; ...
    cin_npref_psth_mn cout_npref_psth_mn; ...
    cin_pref_psth_good_err cout_pref_psth_good_err; ...
    cin_npref_psth_good_err cout_npref_psth_good_err];

pref_combo = [cin_pref_psth_good cout_pref_psth_good];

normsort = @(A,B) sort_by_peak(norm_by_peak(A,B),B);
%normsort = @(A,B) sort_by_peak(A,B);
cax = [0 1]
fh=figure(3); clf

s(1)=subplot(221);
imagesc(normsort(cin_pref_psth_good, pref_combo),'x',cin_t(good_cint));
caxis(cax);
s(2)=subplot(222);
imagesc(normsort(cout_pref_psth_good, pref_combo),'x',cout_t(good_coutt));
caxis(cax);


%colormap(jet)
s(3)=subplot(223);
imagesc(normsort(cin_npref_psth_good, pref_combo),'x',cin_t(good_cint))
caxis(cax);
s(4)=subplot(224);
imagesc(normsort(cout_npref_psth_good, pref_combo),'x',cout_t(good_coutt))
caxis(cax);
%colormap(flipud(gray.^.7));
ylabel(s(3),'cell # (sorted by peak time)')
%ylabel(subplot(223),'cell # (sorted as above)')
hold(s(1),'on')
hold(s(2),'on')
hold(s(3),'on')
hold(s(4),'on')
plot(s(1),[ 0 0], [0 1000],'k')
plot(s(2),[ 0 0], [0 1000],'k')
plot(s(3),[ 0 0], [0 1000],'k')
plot(s(4),[ 0 0], [0 1000],'k')

colormap(colormapLinear([1 1 1].*.0).^.55)

set(subplot(221),'ytick',[300 600])
set(subplot(223),'ytick',[300 600])
set(subplot(222),'ytick',[])
set(subplot(224),'ytick',[])
set(subplot(222),'ytick',[])
set(subplot(224),'ytick',[])
xlabel(subplot(223),'time from stim onset (s)')
xlabel(subplot(224),'from movement (s)')
s(2).YTick = sum(good_cells);
cb = colorbar;
cb.Position = cb.Position + [.1 .35 .02 -.275];
set(cb,'ytick',[0 1],'yticklabel',{'0', 'max'},'fontsize',fsz)
title(cb,'norm FR')
ppos = [8 10 fw 1.5*fht ]
set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]))

for ii = 1:4
    s(ii).TickDir = 'out';
    box(s(ii),'off')
    pbaspect(s(ii),[1 .8 1])
    
end
set(s(1),'XTickLabel',[]); set(s(2),'XTickLabel',[])
%colormap((colormapLinear([1 1 1].*0,50).^.7))
print(fh, fullfile(dp.psth_fig_dir, 'sequence_plot'),...
    '-dsvg', '-painters')


%%
[~,i,j] = sort_by_peak(pref_combo);
pref_combo_t = [cin_t(good_cint) 2+cout_t(good_coutt)];
figure(100); clf
plot(pref_combo_t(j),'.-')
%%
good_cint   = cin_t >= xlim_on(1) & cin_t <= xlim_on(2);
good_coutt  = cout_t >= xlim_off(1) & cout_t <= xlim_off(2);

% good_cint   = cin_t >= -.5 & cin_t <= 1.5;
% good_coutt  = cout_t >= -.5 & cout_t <= .5;

figure(4); clf
ax = subplot(121)
A = cin_psth_hit(good_cells,good_cint,1,2);
B = cout_psth_hit(good_cells,good_coutt,1,2);
C = [A B];
imagesc(normsort(A,C),...
    'x',cin_t(good_cint))
xlim
ax = subplot(122)
imagesc(normsort(B,C),...
    'x',cout_t(good_coutt))

cm = colormapLinear([0 0 0], 49)
colormap(flipud(bone))
colormap(cm)
%colormap(parula)



%%
%%%% PLOT OLD SEQUENCE PLOT
fh1 = figure(1); clf
set(fh1,'paperpositionmode','auto')
subplot(121)
%imagesc(norm_by_peak(psthL(sortA,:)),'x',tA,[0 1]);
imagesc(psth_cin(sortA,:),'x',tA,[0 1])
hold on
plot([0 0], ylim, 'k', 'linewidth', 2);
xlim([-0.25 0.98])
xlabel('Time from stimulus on (s)')
set(gca,'fontsize',12)
ylabel('Cell # (sorted by peak time)')
xticks([0 0.5])
pbaspect([1 1.5 1])
pause(1)

subplot(122)
%imagesc(norm_by_peak(psthR(sortA,:)),'x',tB,[0 1]);
imagesc(psth_cout(sortA,:),'x',tB,[0 1]);
hold on
plot([0 0], ylim, 'k', 'linewidth', 2)
xlim([-0.25 0.98])
xlabel('Time from stimulus off (s)')
pbaspect([1 1.5 1])
xticks([0 0.5])
colormap(flipud(bone.^.4))
%colormap(colormapRedBlue)
%colormap(flipud(colormapBlues.^.5))
set(gca,'fontsize',12)
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 10 10];
fig.PaperPositionMode = 'Manual';
fig.PaperSize = [10 10];

seqplotfname = fullfile(dp.psth_fig_dir,'sequence_plot');
print(fh1, seqplotfname , '-depsc')

corr(sortA, sortB, 'type', 'kendall')

good_cellids = cellids(good_cells);
sorted_cellids = good_cellids(sortA);
prestim = find(sortAt>find(tA>0,1),1);
sorted_stim_cells = sorted_cellids(1:prestim);
sorted_post_stim_cells = sorted_cellids(prestim+1:end);
save(fullfile(dp.ephys_summary_dir, 'sorted_cells.mat'),'sorted_cellids','sorted_stim_cells','sorted_post_stim_cells')
%%

fh = figure(2); clf


ax(1) = subplot(121);hold(ax(1),'on');
ax(2) = subplot(122);hold(ax(2),'on');

plot(ax(1),[ 0 0], [0 100],'k')
plot(ax(2),[ 0 0], [0 100],'k')

xlim(ax(1), [cin_t(find(good_cint,1,'first')) cin_t(find(good_cint,1,'last'))])
xlim(ax(2), [cout_t(find(good_coutt,1,'first')) cout_t(find(good_coutt,1,'last'))])

plot(ax(1),cin_t(good_cint),cin_pref_psth_mn,'color',pref_color,'linewidth',2)
plot(ax(1),cin_t(good_cint),cin_npref_psth_mn,'color',npref_color,'linewidth',2)
plot(ax(1),cin_t(good_cint),cin_pref_psth_good_err,'--','color',pref_color,'linewidth',1)
plot(ax(1),cin_t(good_cint),cin_npref_psth_good_err,'--','color',npref_color,'linewidth',1)


plot(ax(2),cout_t(good_coutt),cout_pref_psth_mn,'color',pref_color,'linewidth',2)
plot(ax(2),cout_t(good_coutt),cout_npref_psth_mn,'color',npref_color,'linewidth',2)
plot(ax(2),cout_t(good_coutt),cout_pref_psth_good_err,'--','color',pref_color,'linewidth',1)
plot(ax(2),cout_t(good_coutt),cout_npref_psth_good_err,'--','color',npref_color,'linewidth',1)

linkaxes(ax,'y')

ylims = ([floor(min(psths(:))*10)/10 ceil(max(psths(:))*10)/10])
%ylims = [.65 1];
ylim(ax(1),ylims)
ylim(ax(2),ylims)


%xlim(ax(1),[-1.1 1.95])
%xlim(ax(2),[-1.95 1.55])

ax(2).YColor = 'w';
ax(1).TickDir = 'out';
ax(2).TickDir = 'out';

%ylabel(ax(1), {'population average' 'normalizeed firing rate'})
ylabel(ax(1), {'normalizeed firing rate'})
xlabel(ax(1), 'time from stim onset (s)')
xlabel(ax(2), 'from movement (s)')
set(fh,'position',[7 5 6 3 ],'papersize',[5 3],'paperpositionmode','auto')

print(fh, fullfile(dp.psth_fig_dir, 'population_psth'),...
    '-dsvg', '-painters')



%%
figure; histogram(prefp)
%%
cm = color_set(2);
fh = figure(10); clf
set(fh,'position',[1 1 3 3], 'papersize', [3 3])
alvl = .05;
pcts = [sum(psr&prefp<alvl) sum(psr==0&prefp<alvl) sum(prefp>=alvl)];
labels = {'right' 'left' 'non-selective'};
p= pie(pcts,labels);
p(1).FaceColor = cm(2,:);
p(1).EdgeColor = [1 1 1];
p(3).FaceColor = cm(1,:);
p(3).EdgeColor = [1 1 1];
p(5).FaceColor = [1 1 1].*.9;
p(5).EdgeColor = [1 1 1];
piechartname = fullfile(dp.psth_fig_dir,'pie_chart');
print(fh, piechartname , '-dsvg','-painters')
