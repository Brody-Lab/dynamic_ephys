close all
clear all
dp = set_dyn_path
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
    if mod(cc,50)==0, fprintf([num2str(cc) '...']); end
    d=dyn_cell_packager(cellids(cc),'datadir',dp.celldat_dir);
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
xlim1 = -.25;
cin_ind = cin_t>-.25 & cin_t<1;%5;
cout_ind = cout_t>xlim1;
save_name = 'sequence_plot';


good_cells = mean(cout_psth(:,cout_ind,1,1),2) > 2;
pref_r = nanmean(cout_psth(:,cout_t > 0,3,1),2) > nanmean(cout_psth(:,cout_t>0,2,1),2);

cout_psth_pref_L_L  = cout_psth(good_cells&~pref_r,:,2,1);
cout_psth_pref_L_R  = cout_psth(good_cells&~pref_r,:,3,1);
cout_psth_pref_R_L  = cout_psth(good_cells&pref_r ,:,2,1);
cout_psth_pref_R_R  = cout_psth(good_cells&pref_r ,:,3,1);

cout_psth_pref_L_LA  = cout_psth(good_cells&~pref_r,:,2,2);
cout_psth_pref_L_RA  = cout_psth(good_cells&~pref_r,:,3,2);
cout_psth_pref_R_LA  = cout_psth(good_cells&pref_r ,:,2,2);
cout_psth_pref_R_RA  = cout_psth(good_cells&pref_r ,:,3,2);

cout_psth_pref_L_LB  = cout_psth(good_cells&~pref_r,:,2,3);
cout_psth_pref_L_RB  = cout_psth(good_cells&~pref_r,:,3,3);
cout_psth_pref_R_LB  = cout_psth(good_cells&pref_r ,:,2,3);
cout_psth_pref_R_RB  = cout_psth(good_cells&pref_r ,:,3,3);

sort_both_choices_together = 1;
if sort_both_choices_together
combo_pref_L = norm_by_peak([cout_psth_pref_L_L(:,cout_ind) cout_psth_pref_L_R(:,cout_ind)]);
combo_pref_R = norm_by_peak([cout_psth_pref_R_L(:,cout_ind) cout_psth_pref_R_R(:,cout_ind)]);
[combo_sortedL, sortL] = sort_by_peak(combo_pref_L);
[combo_sortedR, sortR] = sort_by_peak(combo_pref_R);

combo_pref_LA = norm_by_peak([cout_psth_pref_L_LA(:,cout_ind) cout_psth_pref_L_RA(:,cout_ind)]);
combo_pref_RA = norm_by_peak([cout_psth_pref_R_LA(:,cout_ind) cout_psth_pref_R_RA(:,cout_ind)]);
[combo_sortedLA, sortLA] = sort_by_peak(combo_pref_LA);
[combo_sortedRA, sortRA] = sort_by_peak(combo_pref_RA);

combo_pref_LB = norm_by_peak([cout_psth_pref_L_LB(:,cout_ind) cout_psth_pref_L_RB(:,cout_ind)]);
combo_pref_RB = norm_by_peak([cout_psth_pref_R_LB(:,cout_ind) cout_psth_pref_R_RB(:,cout_ind)]);
[combo_sortedLB, sortLB] = sort_by_peak(combo_pref_LB);
[combo_sortedRB, sortRB] = sort_by_peak(combo_pref_RB);

else
[combo_sortedL, sortL] = sort_by_peak(cout_psth_pref_L_L);
[combo_sortedR, sortR] = sort_by_peak(cout_psth_pref_R_R);
end

corr(sortLA, sortLB, 'type', 'kendall')
corr(sortRA, sortRB, 'type', 'kendall')


psth_L_L  = combo_pref_L(:,1:sum(cout_ind));
psth_L_R  = combo_pref_L(:,sum(cout_ind):end);    
psth_R_L  = combo_pref_R(:,1:sum(cout_ind));
psth_R_R  = combo_pref_R(:,sum(cout_ind):end);

psth_L_LA  = combo_pref_LA(:,1:sum(cout_ind));
psth_L_RA  = combo_pref_LA(:,sum(cout_ind):end);    
psth_R_LA  = combo_pref_RA(:,1:sum(cout_ind));
psth_R_RA  = combo_pref_RA(:,sum(cout_ind):end);

psth_L_LB  = combo_pref_LB(:,1:sum(cout_ind));
psth_L_RB  = combo_pref_LB(:,sum(cout_ind):end);    
psth_R_LB  = combo_pref_RB(:,1:sum(cout_ind));
psth_R_RB  = combo_pref_RB(:,sum(cout_ind):end);


%psth_L_L  = cout_psth_pref_L_L(:,cout_ind);
%psth_L_R  = cout_psth_pref_L_R(:,cout_ind);   
%psth_R_L  = cout_psth_pref_R_L(:,cout_ind);
%psth_R_R  = cout_psth_pref_R_R(:,cout_ind);   

%psth_L_LA  = cout_psth_pref_L_LA(:,cout_ind);
%psth_L_RA  = cout_psth_pref_L_RA(:,cout_ind);   
%psth_R_LA  = cout_psth_pref_R_LA(:,cout_ind);
%psth_R_RA  = cout_psth_pref_R_RA(:,cout_ind);   

%psth_L_LB  = cout_psth_pref_L_LB(:,cout_ind);
%psth_L_RB  = cout_psth_pref_L_RB(:,cout_ind);   
%psth_R_LB  = cout_psth_pref_R_LB(:,cout_ind);
%psth_R_RB  = cout_psth_pref_R_RB(:,cout_ind);   

tB = cout_t(cout_ind);


%%%% PLOT SEQUENCE PLOT
fh1 = figure(1); clf
set(fh1,'paperpositionmode','auto')
subplot(221)
imagesc(psth_L_L(sortL,:),'x',tB,[0 1]);
hold on
plot([0 0], ylim, 'k', 'linewidth', 2);
xlim([xlim1 0.98])
xlabel('Time from stimulus off (s)')
set(gca,'fontsize',12)
ylabel('Cell # (sorted by peak time)')
xticks([0 0.5])
pbaspect([1 1.5 1])
pause(1)

subplot(222)
imagesc(psth_L_R(sortL,:),'x',tB,[0 1]);
hold on
plot([0 0], ylim, 'k', 'linewidth', 2)
xlim([xlim1 0.98])
xlabel('Time from stimulus off (s)')
pbaspect([1 1.5 1])
xticks([0 0.5])
colormap(flipud(bone.^.4))
set(gca,'fontsize',12)
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 10 10];
fig.PaperPositionMode = 'Manual';
fig.PaperSize = [10 10];

subplot(223)
imagesc(psth_R_L(sortR,:),'x',tB,[0 1]);
hold on
plot([0 0], ylim, 'k', 'linewidth', 2);
xlim([xlim1 0.98])
xlabel('Time from stimulus off (s)')
set(gca,'fontsize',12)
ylabel('Cell # (sorted by peak time)')
xticks([0 0.5])
pbaspect([1 1.5 1])
pause(1)

subplot(224)
imagesc(psth_R_R(sortR,:),'x',tB,[0 1]);
hold on
plot([0 0], ylim, 'k', 'linewidth', 2)
xlim([xlim1 0.98])
xlabel('Time from stimulus off (s)')
pbaspect([1 1.5 1])
xticks([0 0.5])
colormap(flipud(bone.^.4))
set(gca,'fontsize',12)
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 10 10];
fig.PaperPositionMode = 'Manual';
fig.PaperSize = [10 10];

print(fullfile(dp.psth_fig_dir, 'choice_sequence'),'-depsc')


%%%% PLOT SEQUENCE PLOT
fh1 = figure(2); clf
set(fh1,'paperpositionmode','auto')
subplot(221)
imagesc(psth_L_LA(sortLA,:),'x',tB,[0 1]);
hold on
plot([0 0], ylim, 'k', 'linewidth', 2);
xlim([xlim1 0.98])
xlabel('Time from stimulus off (s)')
set(gca,'fontsize',12)
ylabel('Cell # (sorted by peak time)')
xticks([0 0.5])
pbaspect([1 1.5 1])
pause(1)

subplot(222)
imagesc(psth_L_RA(sortLA,:),'x',tB,[0 1]);
hold on
plot([0 0], ylim, 'k', 'linewidth', 2)
xlim([xlim1 0.98])
xlabel('Time from stimulus off (s)')
pbaspect([1 1.5 1])
xticks([0 0.5])
colormap(flipud(bone.^.4))
set(gca,'fontsize',12)
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 10 10];
fig.PaperPositionMode = 'Manual';
fig.PaperSize = [10 10];

subplot(223)
imagesc(psth_R_LA(sortRA,:),'x',tB,[0 1]);
hold on
plot([0 0], ylim, 'k', 'linewidth', 2);
xlim([xlim1 0.98])
xlabel('Time from stimulus off (s)')
set(gca,'fontsize',12)
ylabel('Cell # (sorted by peak time)')
xticks([0 0.5])
pbaspect([1 1.5 1])
pause(1)

subplot(224)
imagesc(psth_R_RA(sortRA,:),'x',tB,[0 1]);
hold on
plot([0 0], ylim, 'k', 'linewidth', 2)
xlim([xlim1 0.98])
xlabel('Time from stimulus off (s)')
pbaspect([1 1.5 1])
xticks([0 0.5])
colormap(flipud(bone.^.4))
set(gca,'fontsize',12)
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 10 10];
fig.PaperPositionMode = 'Manual';
fig.PaperSize = [10 10];

print(fullfile(dp.psth_fig_dir, 'choice_sequence_A'),'-depsc')


%%%% PLOT SEQUENCE PLOT
fh1 = figure(3); clf
set(fh1,'paperpositionmode','auto')
subplot(221)
imagesc(psth_L_LB(sortLA,:),'x',tB,[0 1]);
hold on
plot([0 0], ylim, 'k', 'linewidth', 2);
xlim([xlim1 0.98])
xlabel('Time from stimulus off (s)')
set(gca,'fontsize',12)
ylabel('Cell # (sorted by peak time)')
xticks([0 0.5])
pbaspect([1 1.5 1])
pause(1)

subplot(222)
imagesc(psth_L_RB(sortLA,:),'x',tB,[0 1]);
hold on
plot([0 0], ylim, 'k', 'linewidth', 2)
xlim([xlim1 0.98])
xlabel('Time from stimulus off (s)')
pbaspect([1 1.5 1])
xticks([0 0.5])
colormap(flipud(bone.^.4))
set(gca,'fontsize',12)
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 10 10];
fig.PaperPositionMode = 'Manual';
fig.PaperSize = [10 10];

subplot(223)
imagesc(psth_R_LB(sortRA,:),'x',tB,[0 1]);
hold on
plot([0 0], ylim, 'k', 'linewidth', 2);
xlim([xlim1 0.98])
xlabel('Time from stimulus off (s)')
set(gca,'fontsize',12)
ylabel('Cell # (sorted by peak time)')
xticks([0 0.5])
pbaspect([1 1.5 1])
pause(1)

subplot(224)
imagesc(psth_R_RB(sortRA,:),'x',tB,[0 1]);
hold on
plot([0 0], ylim, 'k', 'linewidth', 2)
xlim([xlim1 0.98])
xlabel('Time from stimulus off (s)')
pbaspect([1 1.5 1])
xticks([0 0.5])
colormap(flipud(bone.^.4))
set(gca,'fontsize',12)
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 10 10];
fig.PaperPositionMode = 'Manual';
fig.PaperSize = [10 10];

print(fullfile(dp.psth_fig_dir, 'choice_sequence_B'),'-depsc')

