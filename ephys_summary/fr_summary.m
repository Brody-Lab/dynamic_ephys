close all
clear all

cd ~/Dropbox/spikes/bin/ephys_summary/
addpath ~/ratter/Manuscripts/TimHanks/PBupsPhys/Code/
addpath ~/ratter/Manuscripts/TimHanks/PBupsPhys/Code/Carlosbin
addpath ~/ratter/svn_papers/TimHanks/PBupsPhys/Code/
addpath ~/ratter/svn_papers/TimHanks/PBupsPhys/Code/Carlosbin
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


plot_sorting_data = true;
% decide which timepoints to use for the plot
xlim1 = -.25;
cin_ind = cin_t>0;


avg_fr = mean(cin_psth_hit(:,cin_ind,1,1),2);
good_cells = avg_fr > 2;
avg_fr = avg_fr(good_cells);
mean_population = mean(avg_fr );
mean_population_std = std(avg_fr );

cout_ind = cout_t >= 0 & cout_t <=0.125 ;
avg_fr_choice = mean(cout_psth_hit(good_cells,cout_ind,1,1),2);
mean_population_at_choice = mean(avg_fr_choice);
mean_population_at_choice_std = std(avg_fr_choice);
avg_fr_left_choice = mean(cout_psth_hit(good_cells,cout_ind,2,1),2);
avg_fr_right_choice = mean(cout_psth_hit(good_cells,cout_ind,3,1),2);
avg_lr_diff = avg_fr_right_choice-avg_fr_left_choice;

figure;
hist(avg_fr,50,'k');
ylabel('# Cells')
xlabel('Avg. Firing Rate (Hz)')
set(gca, 'fontsize',16)
pbaspect([1.5 1 1])
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 6];
fig.PaperPositionMode = 'Manual';
fig.PaperSize = [6 6];
print('~/Dropbox/spikes/figures/PSTH/avg_fr','-dsvg')


figure;
hist(avg_lr_diff,50,'k');
ylabel('# Cells')
xlabel('\Delta FR (R - L choice,Hz)')
set(gca, 'fontsize',16)
pbaspect([1.5 1 1])
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 6];
fig.PaperPositionMode = 'Manual';
fig.PaperSize = [6 6];
print('~/Dropbox/spikes/figures/PSTH/avg_fr_diff','-dsvg')


avg_psth_cin = mean(cin_psth(good_cells,:,1,1)); 
avg_psth_cin_L = mean(cin_psth(good_cells,:,2,1)); 
avg_psth_cin_R = mean(cin_psth(good_cells,:,3,1)); 

avg_psth_cout = mean(cout_psth(good_cells,:,1,1)); 
avg_psth_cout_L = mean(cout_psth(good_cells,:,2,1)); 
avg_psth_cout_R = mean(cout_psth(good_cells,:,3,1)); 

figure; 

plot(cin_t,avg_psth_cin,'k','linewidth',2)
hold on
ylim([10 20])
plot([0 0], ylim(), 'k--')
pbaspect([2 1 1])
ylabel('avg fr')
xlabel('time from stim on')
set(gca,'fontsize',16)
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 6];
fig.PaperPositionMode = 'Manual';
fig.PaperSize = [6 6];
print('~/Dropbox/spikes/figures/PSTH/avg_psth_cin','-dsvg')

figure; 
pbaspect([2 1 1])
plot(cout_t,avg_psth_cout,'k','linewidth',2)
ylim([10 20])
hold on;
plot([0 0], ylim(), 'k--')
pbaspect([2 1 1])
ylabel('avg fr')
xlabel('time from stim off')
set(gca,'fontsize',16)

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 6];
fig.PaperPositionMode = 'Manual';
fig.PaperSize = [6 6];
print('~/Dropbox/spikes/figures/PSTH/avg_psth_cout','-dsvg')



