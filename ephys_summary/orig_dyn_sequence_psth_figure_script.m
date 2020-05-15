close all
clear all

addpath ~/ratter/Manuscripts/TimHanks/PBupsPhys/Code/
addpath ~/ratter/Manuscripts/TimHanks/PBupsPhys/Code/Carlosbin
addpath ~/ratter/ExperPort/bin
addpath ~/ratter/Analysis/Pbups
addpath ~/ratter/ExperPort/MySQLUtility
addpath ~/ratter/ExperPort/Analysis
addpath ~/ratter/ExperPort/Analysis/SameDifferent/
addpath ~/ratter/ExperPort/HandleParam
addpath ~/ratter/Analysis/helpers
addpath ~/Dropbox/spikes/cell_packager_data
addpath ~/Dropbox/spikes/bin
addpath ~/Dropbox/spikes/bin/tuning_curves/
%% get full cell_list; *all* single units recorded in pbups project
cell_list = dyn_cells_db;   % That has to be run once to create cell_list
%% compile a bunch of psths aligned to the center poke time and the center out time
% produce a stim on aligned plot and a center out aligned plot
cin_align_ind = 9;
cout_align_ind = 8;

% turn long_trials_only on if you want to only use the long trials for the
% figures
long_trial_dur = 1;
long_trials_only = false;

select_str = 'strcmp(region,''fof'') & normmean > .5' ;
cellids = cell2mat(extracting(cell_list, 'cellid', select_str));
psr = cell2mat(extracting(cell_list, 'prefsideright', select_str));
prefp = cell2mat(extracting(cell_list, 'prefp', select_str));
ncells = size(cellids);

for cc = 1:ncells
    try
    fprintf([num2str(cc) '...'])
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
        cin_psth(cc,:,2,xx)    = nanmean(cin_frates(good&poke_l,:));
        cin_psth(cc,:,3,xx)    = nanmean(cin_frates(good&poke_r,:));
        
        cout_psth(cc,:,1,xx)   = nanmean(cout_frates(good,:));
        cout_psth(cc,:,2,xx)   = nanmean(cout_frates(good&poke_l,:));
        cout_psth(cc,:,3,xx)   = nanmean(cout_frates(good&poke_r,:));
        
        % hits only
        cin_psth_hit(cc,:,1,xx)    = nanmean(cin_frates(hit&good,:));
        cin_psth_hit(cc,:,2,xx)    = nanmean(cin_frates(hit&good&poke_l,:));
        cin_psth_hit(cc,:,3,xx)    = nanmean(cin_frates(hit&good&poke_r,:));
        
        cout_psth_hit(cc,:,1,xx)   = nanmean(cout_frates(hit&good,:));
        cout_psth_hit(cc,:,2,xx)   = nanmean(cout_frates(hit&good&poke_l,:));
        cout_psth_hit(cc,:,3,xx)   = nanmean(cout_frates(hit&good&poke_r,:));
        
        % errors only
        cin_psth_err(cc,:,1,xx)    = nanmean(cin_frates(err&good,:));
        cin_psth_err(cc,:,2,xx)    = nanmean(cin_frates(err&good&poke_l,:));
        cin_psth_err(cc,:,3,xx)    = nanmean(cin_frates(err&good&poke_r,:));
        
        cout_psth_err(cc,:,1,xx)   = nanmean(cout_frates(err&good,:));
        cout_psth_err(cc,:,2,xx)   = nanmean(cout_frates(err&good&poke_l,:));
        cout_psth_err(cc,:,3,xx)   = nanmean(cout_frates(err&good&poke_r,:));
    end
    catch
        disp(cellids(cc))
    end
end

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
save_name = 'sequence_plot';

good_cells = mean(cin_psth_hit(:,cin_ind,1,1),2) > 2;
pref_r = nanmean(cin_psth_hit(:,:,2,1),2) > nanmean(cin_psth_hit(:,:,3,1),2);

cin_psth_hit_pref           = vertcat(cin_psth_hit(good_cells&pref_r,:,2,1), ...
    cin_psth_hit(good_cells&~pref_r,:,3,1));
cout_psth_hit_pref          = vertcat(cout_psth_hit(good_cells&pref_r,:,2,1),...
    cout_psth_hit(good_cells&~pref_r,:,3,1));
cin_psth_hit_pref_xval1     = vertcat(cin_psth_hit(good_cells&pref_r,:,2,2), ...
    cin_psth_hit(good_cells&~pref_r,:,3,2));
cout_psth_hit_pref_xval1    = vertcat(cout_psth_hit(good_cells&pref_r,:,2,2), ...
    cout_psth_hit(good_cells&~pref_r,:,3,2));
cin_psth_hit_pref_xval2     = vertcat(cin_psth_hit(good_cells&pref_r,:,2,3),...
    cin_psth_hit(good_cells&~pref_r,:,3,3));
cout_psth_hit_pref_xval2    = vertcat(cout_psth_hit(good_cells&pref_r,:,2,3),...
    cout_psth_hit(good_cells&~pref_r,:,3,3));

combo_psthA = [cin_psth_hit_pref_xval1(:,cin_ind) cout_psth_hit_pref_xval1(:,cout_ind)];
combo_psthB = [cin_psth_hit_pref_xval2(:,cin_ind) cout_psth_hit_pref_xval2(:,cout_ind)];
%combo_psthA = [cin_psth_hit_pref_xval1(:,:) cout_psth_hit_pref_xval1(:,:)];
%combo_psthB = [cin_psth_hit_pref_xval2(:,:) cout_psth_hit_pref_xval2(:,:)];

[combo_sortedA, sortA] = sort_by_peak(combo_psthA);
[combo_sortedB, sortB] = sort_by_peak(combo_psthB);
if plot_sorting_data
    psth_all = [cin_psth_hit_pref_xval1(:,cin_ind) cout_psth_hit_pref_xval1(:,cout_ind)];
%    psth_all = [cin_psth_hit_pref_xval1(:,:) cout_psth_hit_pref_xval1(:,:)];
    for i=1:size(psth_all,1)
        psth_all(i,:) = norm_by_peak(psth_all(i,:));
    end 

    psthL = psth_all(:,1:sum(cin_ind));
    psthR = psth_all(:,sum(cin_ind)+1:end);
%    psthL = psth_all(:,1:120);
%    psthR = psth_all(:,121:end);
%    psthL = psthL(:,cin_ind);
%    psthR = psthR(:,cout_ind);    

else
    psthL = cin_psth_hit_pref_xval2(:,cin_ind);
    psthR = cout_psth_hit_pref_xval2(:,cout_ind);
end
    
tA = cin_t(cin_ind);
tB = cout_t(cout_ind);


%%%% PLOT SEQUENCE PLOT
fh1 = figure(1); clf
set(fh1,'paperpositionmode','auto')
subplot(121)
%imagesc(norm_by_peak(psthL(sortA,:)),'x',tA,[0 1]);
imagesc(psthL(sortA,:),'x',tA,[0 1])
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
imagesc(psthR(sortA,:),'x',tB,[0 1]);
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

print(fh1,'~/Dropbox/spikes/figures/PSTH/sequence_plot' , '-depsc')

corr(sortA, sortB, 'type', 'kendall')

good_cellids = cellids(good_cells);
sorted_cellids = good_cellids(sortA);
sorted_stim_cells = sorted_cellids(1:60);
sorted_post_stim_cells = sorted_cellids(61:end);
save('~/Dropbox/spikes/bin/ephys_summary/sorted_cells.mat','sorted_cellids','sorted_stim_cells','sorted_post_stim_cells')


keyboard
%% PLOT AVERAGE PSTH
% redo this so that the cells actually line up with preferred non-preferred


pref_r = nanmean(cin_psth_hit(:,:,2,1),2) > nanmean(cin_psth_hit(:,:,3,1),2);
switch 0
    case 0
        cell_to_plot = 16857;
%        cell_to_plot = 17784;
%        cell_to_plot = 18181;
        good_cells = cellids == cell_to_plot;
        titlestr = ['cell ' num2str(cell_to_plot)];
    case 1
        good_cells = nanmean(cin_psth(:,:,1,1),2) > 10;
        titlestr = 'all cells with FR > 10';
end
cin_psth_hit_pref = nanmean(vertcat(cin_psth_hit(good_cells&pref_r,:,2,1), cin_psth_hit(good_cells&~pref_r,:,3,1)),1);
cin_psth_err_pref = nanmean(vertcat(cin_psth_err(good_cells&pref_r,:,2,1), cin_psth_err(good_cells&~pref_r,:,3,1)),1);
cin_psth_hit_nonpref = nanmean(vertcat(cin_psth_hit(good_cells&~pref_r,:,2,1), cin_psth_hit(good_cells&pref_r,:,3,1)),1);
cin_psth_err_nonpref = nanmean(vertcat(cin_psth_err(good_cells&~pref_r,:,2,1), cin_psth_err(good_cells&pref_r,:,3,1)),1);
cout_psth_hit_pref = nanmean(vertcat(cout_psth_hit(good_cells&pref_r,:,2,1), cout_psth_hit(good_cells&~pref_r,:,3,1)),1);
cout_psth_err_pref = nanmean(vertcat(cout_psth_err(good_cells&pref_r,:,2,1), cout_psth_err(good_cells&~pref_r,:,3,1)),1);
cout_psth_hit_nonpref = nanmean(vertcat(cout_psth_hit(good_cells&~pref_r,:,2,1), cout_psth_hit(good_cells&pref_r,:,3,1)),1);
cout_psth_err_nonpref = nanmean(vertcat(cout_psth_err(good_cells&~pref_r,:,2,1), cout_psth_err(good_cells&pref_r,:,3,1)),1);

% figure out color
cmap = color_set(2);
if pref_r(good_cells)
    nonpref_col = cmap(1,:);
    pref_col = cmap(2,:);
else
    nonpref_col = cmap(2,:);
    pref_col = cmap(1,:);
end


ylims = [0 1.1*max([max(cin_psth_hit_pref); max(cin_psth_hit_nonpref); max(cin_psth_err_pref); max(cin_psth_err_nonpref);max(cout_psth_hit_pref); max(cout_psth_hit_nonpref); max(cout_psth_err_pref); max(cout_psth_err_nonpref)])];
fh2 = figure(2); clf
subplot(1,2,1); hold on
plot(cin_t,cin_psth_hit_pref,   '-' ,'linewidth',2,'color',pref_col)
plot(cin_t,cin_psth_hit_nonpref,'-' ,'linewidth',2,'color',nonpref_col)
plot(cin_t,cin_psth_err_pref,   '--','linewidth',2,'color',pref_col)
plot(cin_t,cin_psth_err_nonpref,'--','linewidth',2,'color',nonpref_col)
xlabel('Time from stim on (s)')
ylabel('mean FR (Hz)')
set(gca,'fontsize',16)
pbaspect([1.5 1 1])
xlim([-1 2])
ylim(ylims)

subplot(1,2,2); hold on;
plot(cout_t,cout_psth_hit_pref,     '-' ,'linewidth',2,'color',pref_col)
plot(cout_t,cout_psth_hit_nonpref,  '-' ,'linewidth',2,'color',nonpref_col)
plot(cout_t,cout_psth_err_pref,     '--','linewidth',2,'color',pref_col)
plot(cout_t,cout_psth_err_nonpref,  '--','linewidth',2,'color',nonpref_col)
xlabel('Time from stim off (s)')
ylabel('mean FR (Hz)')
set(gca,'fontsize',16)
xlim([-2.5 1])
ylim(ylims)
if pref_r(good_cells)
    hl=legend('go R hit','go L hit',' go R error','go L error','location','northwest')
else
    hl=legend('go L hit','go R hit',' go L error','go R error','location','northwest')
end
box(hl,'off')
plot(subplot(121),[0 0], ylim, 'k--','linewidth',1)
plot(subplot(122),[0 0], ylim, 'k--','linewidth',1)
pbaspect([1.5 1 1])
%suptitle(titlestr)
set(fh2, 'paperorientation','landscape')
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 10 10];
fig.PaperPositionMode = 'Manual';
fig.PaperSize = [10 10];
print(fh2, ['~/Dropbox/spikes/figures/PSTH/cell_' num2str(cell_to_plot) ], '-dsvg')

