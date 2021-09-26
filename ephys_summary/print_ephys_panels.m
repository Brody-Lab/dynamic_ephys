close all
clear all
dp = set_dyn_path(1)
[align_strs, align_args]    = dyn_align_LUT;
cout_align_ind  = find(ismember(align_strs,'cpokeout'));
cin_align_ind   = find(ismember(align_strs,'stimstart-cout-mask'));
assert(all([~isempty(cout_align_ind),~isempty(cin_align_ind)]));

% turn long_trials_only on if you want to only use the long trials for the
% figures
long_trial_dur = 1;
long_trials_only = false;
ploterrorbar = 1;
coutstr = 'cpokeout';
repack  = 0;
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
        d=dyn_cell_packager(cellids(cc),'repack',repack);
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
fprintf('\nnumber of cells: %i \nnumber of sessions: %i\n',...
    length(cellids), length(unique(sessids)));

rats = dp.ratlist';

ncells_per_rat = nan(size(rats));
nsess = nan(size(rats));
for rr = 1:length(rats)
    idx = sum(ratnames == rats{rr},2)==4;
    ncells_per_rat(rr) = sum(idx);
    nsess(rr) = length(unique(sessids(idx)));
end
fprintf('\nmin/max cells: %i/%i \nmin/max sessions: %i/%i\n',...
    min(ncells_per_rat),max(ncells_per_rat),min(nsess),max(nsess))
fprintf('\nmean cells: %.2f \nmean sessions: %.2f\n',...
    mean(ncells_per_rat),mean(nsess))
table(rats,ncells_per_rat,nsess)
ncells = sum(ncells_per_rat);

cin_ind         = cin_t > -.25 & cin_t < 2;
good_cells1     = mean(cin_psth_hit(:,cin_ind,1,1),2) > 2;
mn_fr           = nanmean(cin_psth(:, cin_t > 0, 1, 1),2);
active_cells    = mn_fr > normmnth;
sig_cells       = prefp < prefpth;
good_cells      = active_cells & sig_cells;
fprintf(['n good cells = ' num2str(sum(good_cells))])
fprintf('\n%.1f %% (%i/%i) of active cells are signficant\n',...
    100*sum(good_cells)/sum(active_cells),sum(sig_cells&active_cells),sum(active_cells))

%% plot example cells in fig 2B
ppos = [8 10 fw fht ];
cellid = 18181;
[fh, ax] = example_cell_psth('cells',cellid,...
    'cintrange',xlim_on,'couttrange',xlim_off,...
    'coutstr',coutstr,'fig_num',2, ...
    'ploterrorbar',ploterrorbar,'repack',repack,'errorbarfun',@nansem);

ylim(ax,[14 60])
xlim(ax(1),xlim_on)
xlim(ax(2),xlim_off)
ax(1).YTick = [15:15:60];
ax(2).YTick = [15:15:60];
ax(2).YTickLabel = {};

set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],...
    'papersize',ppos([3 4]))
pbaspect(ax(1),[1 .8 1])
pbaspect(ax(2),[1 .8 1])

print(fh, fullfile(dp.psth_fig_dir, ['cell_' num2str(cellid) ]),...
    '-dsvg', '-painters')
%%
cellid = 16857;
[fh, ax] = example_cell_psth('cells',cellid,...
    'cintrange',xlim_on,'couttrange',xlim_off,...
    'coutstr',coutstr,'fig_num',2, ...
    'ploterrorbar',ploterrorbar,'repack',repack,'errorbarfun',@nansem);
ylim(ax,[0 22.5])
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
[fh, ax] = example_cell_psth('cells',cellid,...
    'cintrange',xlim_on,'couttrange',xlim_off,...
    'coutstr',coutstr,'fig_num',2, ...
    'ploterrorbar',ploterrorbar,'repack',repack,'errorbarfun',@nansem);
ylim(ax,[0 42.5])
xlim(ax(1),xlim_on)
xlim(ax(2),xlim_off)
ax(2).YTickLabel = {};
pbaspect(ax(1),[1 .8 1])
pbaspect(ax(2),[1 .8 1])
set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]))
print(fh, fullfile(dp.psth_fig_dir, ['cell_' num2str(cellid) ]),...
    '-dsvg', '-painters')

%% raster plots for supplement
ex_cellids = [18181 16857 17784];
max_nt = 200;
for cc = 1:length(ex_cellids)
    ppos_ = [8 10 1.75*fw 1.25*fht ];
    cellid = ex_cellids(cc);
    [fh, ax] = example_cell_raster('cells',cellid,...
        'cintrange',xlim_on,'couttrange',xlim_off,...
        'coutstr',coutstr,'fig_num',2,'repack',repack,'max_ntrials',max_nt);
    set(ax(1),'ZLim',get(ax(2),'Zlim'));
    set(fh,'position',ppos_,'paperposition',[0 0 ppos_([3 4])],...
        'papersize',ppos_([3 4]));
    print(fh, fullfile(dp.psth_fig_dir, ['cell_' num2str(cellid) '_raster' ]),...
        '-dsvg', '-painters');
end

%% plot side x outcome population psth for supplement
[fh, ax] = example_cell_psth('cells',cellids(good_cells),...
    'ploterrorbar',1,'meta',1,...
    'cintrange',xlim_on,'couttrange',xlim_off,...
    'coutstr',coutstr,'fig_num',2, 'norm', 'onset', 'repack',repack)
xlim(ax(1),xlim_on)
xlim(ax(2),xlim_off)

ax(2).YTickLabel = {};
pbaspect(ax(1),[1 .8 1])
pbaspect(ax(2),[1 .8 1])
set(fh,'position',ppos,'paperposition',[0 0 ppos([3 4])],'papersize',ppos([3 4]))
print(fh, fullfile(dp.psth_fig_dir, 'selective_cells_psth'),...
    '-dsvg', '-painters')

%% plot chronometric population psth for fig 2D
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
    %%
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

%% plot sorted auc fig 2C
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
    ylabel({'Cell # (sorted by latency)'})

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

