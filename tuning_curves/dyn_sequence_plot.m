addpath('~/ratter/svn_papers/TimHanks/PBupsPhys/Code/');
addpath('~/ratter/svn_papers/TimHanks/PBupsPhys/Code/Carlosbin');
%% get full cell_list; *all* single units recorded in pbups project
cell_list = dyn_cells_db('force',0);   % That has to be run once to create cell_list
%%
align_ind = 14;
if align_ind == 8
    xlab = 'time from center out (s)';
elseif align_ind == 9
    xlab = 'time from stim on (s)';
elseif align_ind == 14
    xlab = 'time from stim on (s)';
else
    xlab = '';
end
select_str = 'strcmp(region,''fof'')' ;
cellids = cell2mat(extracting(cell_list, 'cellid', select_str));
psr = cell2mat(extracting(cell_list, 'prefsideright', select_str));
prefp = cell2mat(extracting(cell_list, 'prefp', select_str));
ncells = size(cell_list,1)-1;
for cc = 1:ncells
    fprintf([num2str(cc) '...'])
    d=dyn_cell_packager(cellids(cc));
    this_frates = d.frate{align_ind};
      
    
    if cc == 1
        t = d.frate_t{align_ind};
        bin_size = diff(t([1 2]));
        ntp = length(t);
        pref_psth = nan(length(cellids),ntp);
        nonpref_psth = nan(length(cellids),ntp);
        pref_psth_hit = nan(length(cellids),ntp);
        nonpref_psth_hit = nan(length(cellids),ntp);
        pref_psth_err = nan(length(cellids),ntp);
        nonpref_psth_err = nan(length(cellids),ntp);
    end
    
    if 0 
        figure(2);
        subplot(211)
        imagesc(this_frates,'x',t); caxis([-1 10]); colormap([.5 0 0; bone]);
        subplot(212)
        plot(t,nanmean(this_frates))
        pause(.5)
    end
    
    poke_r = d.trials.rat_dir==1;
    poke_l = d.trials.rat_dir==-1;
    hit = d.trials.hit==1;
    stim_dur = d.trials.cpoke_end - d.trials.stim_start;
    good = stim_dur > 1.2;
    
    if psr(cc) == 1
        pref_psth_hit(cc,:) = nanmean(this_frates(good&poke_r&hit,:));
        nonpref_psth_hit(cc,:) = nanmean(this_frates(good&poke_l&hit,:));
        pref_psth_err(cc,:) = nanmean(this_frates(good&poke_r&hit==0,:));
        nonpref_psth_err(cc,:) = nanmean(this_frates(good&poke_l&hit==0,:));
    else
        pref_psth_hit(cc,:) = nanmean(this_frates(good&poke_l&hit,:));
        nonpref_psth_hit(cc,:) = nanmean(this_frates(good&poke_r&hit,:));
        pref_psth_err(cc,:) = nanmean(this_frates(good&poke_l&hit==0,:));
        nonpref_psth_err(cc,:) = nanmean(this_frates(good&poke_r&hit==0,:));
    end
    
end
%%

figure(1); clf
good_cells = true(size(pref_psth_hit(:,1)));
good_cells = nanmean(pref_psth_hit,2)>2 & prefp < .05;
pref_psth_trim = pref_psth_hit(good_cells,:);
nonpref_psth_trim = nonpref_psth_hit(good_cells,:);
pref_psth_trim_err = pref_psth_err(good_cells,:);
nonpref_psth_trim_err = nonpref_psth_err(good_cells,:);

[sorted_pref_psth, sort_order] = sort_by_peak(pref_psth_trim);
sorted_pref_psth = pref_psth_trim(sort_order,:);

subplot(231);
imagesc((norm_by_peak(pref_psth_trim(sort_order,:))),'x',t,[0 1])
hold on
plot([0 0],ylim,'k')    

title('pref hits')
xlabel(xlab)
ylabel({'cell #' '(sorted by peak time in preferred trials'})
set(gca,'linewidth',1.5,'fontsize',15)
colorbar

subplot(232);
imagesc(norm_by_peak(nonpref_psth_trim(sort_order,:),sorted_pref_psth),'x',t,[0 1])
hold on
plot([0 0],ylim,'k')

colormap(parula)
title('non-pref hits')
xlabel(xlab)
ylabel({'cell #' '(sorted by peak time in preferred trials'})
set(gca,'linewidth',1.5,'fontsize',15)
colorbar

subplot(234);
imagesc(norm_by_peak(pref_psth_trim_err(sort_order,:),sorted_pref_psth),'x',t,[0 1])
hold on
plot([0 0],ylim,'k')

title('went pref in error')
xlabel(xlab)
ylabel({'cell #' '(sorted by peak time in preferred trials'})
set(gca,'linewidth',1.5,'fontsize',15)
colorbar

subplot(235);
imagesc(norm_by_peak(nonpref_psth_trim_err(sort_order,:),sorted_pref_psth),'x',t,[0 1])
hold on
plot([0 0],ylim,'k')

colormap(parula)
title('went non-pref in error')
xlabel(xlab)
ylabel({'cell #' '(sorted by peak time in preferred trials'})
%colormap(flipud(colormapBlues))

set(gca,'linewidth',1.5,'fontsize',15)
colorbar

subplot(233);
hit_diff = (norm_by_peak(pref_psth_trim(sort_order,:)))-...
    (norm_by_peak(nonpref_psth_trim(sort_order,:),sorted_pref_psth));
err_diff = (norm_by_peak(pref_psth_trim_err(sort_order,:),sorted_pref_psth))-...
    (norm_by_peak(nonpref_psth_trim_err(sort_order,:),sorted_pref_psth));
imagesc(hit_diff,'x',t,[-1 1])
hold on
plot([0 0],ylim,'k')
title('pref-nonpref (hits)')
xlabel(xlab)
ylabel({'cell #' '(sorted by peak time in preferred trials'})
set(gca,'linewidth',1.5,'fontsize',15)
colorbar

subplot(236);
imagesc(err_diff,'x',t,[-1 1])
hold on
plot([0 0],ylim,'k')
colormap(parula)
title('pref-nonpref (errors)')
xlabel(xlab)
ylabel({'cell #' '(sorted by peak time in preferred trials'})
%colormap(flipud(colormapBlues))
colorbar
colormap(colormapRedBlue.^.6)
%colormap(flipud(bone))
set(gca,'linewidth',1.5,'fontsize',15)

%%
figure(2); clf
r_col = [.75 .2 .2];
g_col = [.2 .75 .2];
ind = sort_order;
plot(t,nanmean(pref_psth_trim(ind,:),1),'-','linewidth',2,'color',g_col)
hold on
plot(t,nanmean(pref_psth_trim_err(ind,:),1),'-','linewidth',2,'color',r_col)

plot(t,nanmean(nonpref_psth_trim(ind,:),1),'--','linewidth',2,'color',g_col)
hold on
plot(t,nanmean(nonpref_psth_trim_err(ind,:),1),'--','linewidth',2,'color',r_col)

legend('pref hit','went pref in error','non-pref hit','went non-pref in error',...
    'location','northwest')
plot([0 0], ylim, 'k','linewidth',1.5)

xlabel(xlab)
ylabel({'mean FR across cells with pre-movement selectivity' 'and <firing rate> > 2'})
%ylabel({'mean FR across cells with' '<firing rate> > 2'})
set(gca,'linewidth',1.5,'fontsize',17)
%%

[prefp_ordered, prefp_order] = sort(prefp);
for cc = 1:length(cellids)
    figure(3); clf
    d = dyn_cell_packager(cellids(prefp_order(cc)));
    subplot(121)
    good = d.trials.rat_dir==d.trials.correct_dir;
    stim_fr_t = d.frate_t{9};
    stim_fr = d.frate{9};
    hold on; plot(stim_fr_t,nanmean(stim_fr(good&d.trials.bup_diff>30,:)),'color',[.5 0 0],'linewidth',2)
    plot(stim_fr_t,nanmean(stim_fr(good&d.trials.bup_diff<30&d.trials.bup_diff>0,:)),...
        'color',[.75 .4 .4 ],'linewidth',2)
    plot(stim_fr_t,nanmean(stim_fr(good&d.trials.bup_diff>-30&d.trials.bup_diff<0,:)),...
        'color',[.4 .4 .75],'linewidth',2)
    plot(stim_fr_t,nanmean(stim_fr(good&d.trials.bup_diff<-30,:)),...
        'color',[0 0 .5],'linewidth',2)
    plot([ 0 0],ylim,'k')
    xlabel('time from stim on (s)')
    ylabel('spks/s')
    
    subplot(122)
    hold on; plot(d.frate_t{8},nanmean(d.frate{8}(good&d.trials.bup_diff>30,:)),'color',[.5 0 0],'linewidth',2)
    plot(d.frate_t{8},nanmean(d.frate{8}(good&d.trials.bup_diff<30&d.trials.bup_diff>0,:)),'color',[.75 .4 .4 ],'linewidth',2)
    plot(d.frate_t{8},nanmean(d.frate{8}(good&d.trials.bup_diff>-30&d.trials.bup_diff<0,:)),'color',[.4 .4 .75],'linewidth',2)
    plot(d.frate_t{8},nanmean(d.frate{8}(good&d.trials.bup_diff<-30,:)),'color',[0 0 .5],'linewidth',2)
    legend({'easy R', 'hard R', 'hard L', 'easy L'})
    xlabel('time from center out (s)')
    ylabel('spks/s')
    title([num2str(cellids(prefp_order(cc))) ' prefp = ' num2str(prefp_ordered(cc))]);
    plot([ 0 0],ylim,'k')
    
    linkaxes([subplot(121) subplot(122)],'y')
    pause
end
%%
figure;
plot(norm_by_peak(pref_psth_trim(1:10,:))')

%% assemble l and r psth
regions = {'fof'};
[align_strs, align_args] = dyn_align_LUT;
align_str = align_strs{align_ind};
save_plots = 0;
minmaxfr = 0;
save_mat = 0;
for rr = 1:length(regions)

    
    select_str = ['strcmp(region,''' regions{rr} ''')  & lr>.75 & mingamma<1.5 '];
    cellids = cell2mat(extracting(cell_list, 'cellid', select_str));
    psr = cell2mat(extracting(cell_list, 'prefsideright', select_str));
    ratnames = cell2mat(extracting(cell_list, 'prefsideright', select_str));
    sessid =  extracting(cell_list,'sessid',select_str);
    prefp =  extracting(cell_list,'prefp',select_str);
    fprintf(['working on cell...' repmat(' ', 1, 14)])
    for cc = 1:length(cellids)
        fprintf([repmat('\b', 1, 14) '%3i ... of %3i'],cc,length(cellids))
        d = cell_packager(cellids(cc));
        this_frates = d.frate{align_ind};
        if cc == 1
            t = d.frate_t{align_ind};
            bin_size = diff(t([1 2]));
            ntp = length(t);
            pref_psth = nan(length(cellids),ntp);
            nonpref_psth = nan(length(cellids),ntp);
            auc   = nan(length(cellids),ntp);
            auc_p = nan(length(cellids),ntp);
        end
        
        poke_r = d.trials.rat_dir==1;
        poke_l = d.trials.rat_dir==-1;
        
        if psr(cc) == 1
            pref_psth(cc,:) = nanmean(this_frates(poke_r,:));
            nonpref_psth(cc,:) = nanmean(this_frates(poke_l,:));
        else
            pref_psth(cc,:) = nanmean(this_frates(poke_l,:));
            nonpref_psth(cc,:) = nanmean(this_frates(poke_r,:));
        end
        %tic
        for tt = 1:length(t)
            [auc(cc,tt), auc_p(cc,tt), this_ci] = ...
                bootroc(this_frates(poke_r,tt),this_frates(poke_l,tt), 500, 95);
        end
        % toc
        
    end
    %%
%     rr = 1;
%     load(fullfile('~/projects/ppcfof_tuning/sequences/',...
%     ['all_' regions{rr} '_aucs_psths_align_' align_str '.mat']));
        f = figure(rr);
    clf(f)
    ind = max(pref_psth,[],2) > minmaxfr;
    
    pref_psth_norm = pref_psth(ind,:) ./ max(pref_psth(ind,:),[],2);
    nonpref_psth_norm = nonpref_psth(ind,:) ./ max(pref_psth(ind,:),[],2);
    
    auc_subset = auc(ind,:);
    auc_p_subset = auc_p(ind,:);
    
    cumsig = cumsum(abs(auc_p_subset-.5)>.499,2)==ceil(.15/bin_size) ;
    cumsig(:,end)=ceil(.15/bin_size);
    [auc_sort_i, auc_sort_j] = find(cumsum(cumsum(cumsig,2),2)==1);
    
    %[~, auc_sort_i] = sort_by_peak(auc_subset);
    
    [pref_norm_sorted, psth_sort_i, psth_sort_j] = sort_by_peak(pref_psth_norm);
    nonpref_norm_sorted = nonpref_psth(psth_sort_i,:) ...
        ./ max(pref_psth(psth_sort_i,:),[],2);
    
    suptitle([select_str ' & maxfr >' num2str(minmaxfr) 'Hz'])
    
    subplot(231,'parent',f)
    imagesc(abs(auc_subset(auc_sort_i,:)-.5),'x',t,[0 .5])
    hold on
    plot([0 0],ylim,'w','linewidth',1.5)
    xlabel(align_str)
    ylabel('cell # (sorted by selectivity onset)')
    title('selectivity')
    colorbar
    
    subplot(232,'parent',f)
    imagesc(pref_psth_norm(auc_sort_i,:), 'x', t,[0 1])
    hold on
    plot([0 0],ylim,'w','linewidth',1.5)
    xlabel(align_str)
    ylabel('cell # (sorted by selectivity onset)')
    title('preferred side psth')
    colorbar
    
    subplot(233,'parent',f)
    imagesc(nonpref_psth_norm(auc_sort_i,:), 'x', t,[0 1]);
    hold on
    plot([0 0],ylim,'w','linewidth',1.5)
    title('non-preferred side psth')
    xlabel(align_str)
    ylabel('cell # (sorted by selectivity onset)')
    colorbar
    
    subplot(234,'parent',f)
    imagesc(abs(auc_subset(psth_sort_i,:)-.5),'x',t, [0 .5])
    hold on
    plot([0 0],ylim,'w','linewidth',1.5)
    xlabel(align_str)
    ylabel('cell #  (sorted by psth peak time)')
    title('selectivity')
    colorbar
    
    subplot(235,'parent',f)
    imagesc(pref_norm_sorted, 'x', t, [0 1])
    hold on
    plot([0 0],ylim,'w','linewidth',1.5)
    xlabel(align_str)
    ylabel('cell # (sorted by psth peak time)')
    title('preferred side psth')
    colorbar
    
    subplot(236,'parent',f)
    imagesc(nonpref_norm_sorted, 'x', t, [0 1]);
    hold on
    plot([0 0],ylim,'w','linewidth',1.5)
    xlabel(align_str)
    ylabel('cell # (sorted by psth peak time)')
    title('non-preferred side psth')
    colorbar
    colormap(hot.^.75)
    
    %
    
    
    %     print(f, fullfile('~/projects/ppcfof_tuning/sequences/',...
    %         [regions{rr} '_premvmt_seq_align_' align_str '.jpg']),'-djpeg');
    %
    %     save(fullfile('~/projects/ppcfof_tuning/sequences/',...
    %         [regions{rr} '_premvmt_seq_align_' align_str '.mat']), 'auc','auc_p','pref_psth',...
    %         'nonpref_psth','cellids','select_str','psr', 'select_str','align_str',...
    %         't','ratnames','sessid','n_trials')
    
    %
    if save_plot
    print(f, fullfile('~/projects/fof_dyn',...
        ['dyn_all_' regions{rr} '_aucs_psths_seq_align_' align_str '.jpg']),'-djpeg');
    end
    if save_mat
    save(fullfile('~/projects/fof_dyn',...
        ['dyn_all_' regions{rr} '_aucs_psths_align_' align_str '.mat']), 'auc','auc_p','pref_psth',...
        'nonpref_psth','cellids','select_str','psr', 'select_str','align_str',...
        't','ratnames','sessid','n_trials')
    end
    % want to know if cell
    
end

%% load and plot sequence for a given
rr = 3;
f = figure(1);
clf(f)
load(    fullfile('~/projects/ppcfof_tuning/sequences/',...
    [regions{rr} '_premvmt_seq_align_' align_str '.mat']), 'auc','auc_p','pref_psth',...
    'nonpref_psth','cellids','select_str','psr','select_str','align_str',...
    't','ratnames','sessid','n_trials')

clf(f)

fr_ind = max(pref_psth,[],2) > minmaxfr;
[~, psth_sort_i, psth_sort_j] = sort_by_peak(pref_psth);

[~, maxt_ind] = max(pref_psth,[],2);

sort_inds = sub2ind(size(pref_psth), psth_sort_i, psth_sort_j);

% can't directly use psth_sort_j for subset plot, because sorted
% differently, in fact, should be worried about all this sorting

ind = abs(auc_p(sort_inds)-.5) > .4999 & fr_ind &  t(maxt_ind)' < 0;
% ;%& fr_ind &;

%[sig_pk_ind,plot_j] = ind2sub(size(pref_psth_norm), sort_inds(sig_peaks));

pref_psth_norm = pref_psth(ind,:) ./ max(pref_psth(ind,:),[],2);
nonpref_psth_norm = nonpref_psth(ind,:) ./ max(pref_psth(ind,:),[],2);
auc_subset = auc(ind,:);
auc_p_subset = auc_p(ind,:);

[pref_norm_sorted, psth_sort_i, psth_sort_j] = sort_by_peak(pref_psth_norm);
nonpref_norm_sorted = nonpref_psth(psth_sort_i,:) ...
    ./ max(pref_psth(psth_sort_i,:),[],2);

cumsig = cumsum(abs(auc_p_subset-.5)>.499,2)==ceil(.15/bin_size) ;
cumsig(:,end)=ceil(.15/bin_size);
[auc_sort_i, auc_sort_j] = find(cumsum(cumsum(cumsig,2),2)==1);

suptitle([select_str ' & maxfr >' num2str(minmaxfr) 'Hz'])

subplot(231,'parent',f)
imagesc(abs(auc_subset(auc_sort_i,:)-.5),'x',t,[0 .5])
hold on
plot([0 0],ylim,'w','linewidth',1.5)
xlabel('time from cout (s)')
ylabel('cell # (sorted by selectivity onset)')
title('selectivity')
colorbar

subplot(232,'parent',f)
imagesc(pref_psth_norm(auc_sort_i,:), 'x', t,[0 1])
hold on
plot([0 0],ylim,'w','linewidth',1.5)
xlabel('time from cout (s)')
ylabel('cell # (sorted by selectivity onset)')
title('preferred side psth')
colorbar

subplot(233,'parent',f)
imagesc(nonpref_psth_norm(auc_sort_i,:), 'x', t,[0 1]);
hold on
plot([0 0],ylim,'w','linewidth',1.5)
title('non-preferred side psth')
xlabel('time from cout (s)')
ylabel('cell # (sorted by selectivity onset)')
colorbar

subplot(234,'parent',f)
imagesc(abs(auc_subset(psth_sort_i,:)-.5),'x',t, [0 .5])
hold on
plot([0 0],ylim,'w','linewidth',1.5)
xlabel('time from cout (s)')
ylabel('cell #  (sorted by psth peak time)')
title('selectivity')
colorbar

subplot(235,'parent',f)
imagesc(pref_norm_sorted, 'x', t, [0 1])
hold on
plot([0 0],ylim,'w','linewidth',1.5)
xlabel('time from cout (s)')
ylabel('cell # (sorted by psth peak time)')
title('preferred side psth')
colorbar

subplot(236,'parent',f)
imagesc(nonpref_norm_sorted, 'x', t, [0 1]);
hold on
plot([0 0],ylim,'w','linewidth',1.5)
xlabel('time from cout (s)')
ylabel('cell # (sorted by psth peak time)')
title('non-preferred side psth')
colorbar
colormap(hot.^.75)


print(f, fullfile('~/projects/ppcfof_tuning/sequences/',...
    [regions{rr} '_premvmt_sigpeak_stim_seq.jpg']),'-djpeg');





%%
[sim_raster, sim_rate, t_filt] = simulate_spikes_rates(avg_pref_fr, t, sum(n_trials));
figure(3); clf
imagesc(sim_rate)
%%
ind = abs(auc_p(sort_inds)-.5) > .4999 & fr_ind &  t(maxt_ind)' < 0;
%ind = abs(auc_p(sort_inds)-.5) > .4999 & fr_ind ;
real_psth = pref_psth(ind,:);
n_cells = size(real_psth,1);
this_n = round(n_trials(ind)/10);
avg_pref_fr = nanmean(real_psth);
figure(12);clf

[sim_raster_full, sim_rate_full, t_filt] = simulate_spikes_rates(avg_pref_fr, t, sum(this_n));
noise_sd = 0;
sim_psth = nan(n_cells,length(t_filt));
sim_rate_assembled = nan(size(sim_rate_full));
%this_fr = nan(n_cells,length(t));
noisesd = .5;
noisesd = .01;
for cc = 1:n_cells
    noisy_fr = avg_pref_fr .* ( 1 + noisesd .* randn(this_n(cc),length(t)));
    %noisy_fr = avg_pref_fr .* ( 1 + noisesd .* randn(1,length(t)));
    noisy_scaled_fr = noisy_fr .* max(real_psth(cc,:))./max(avg_pref_fr);
    this_fr = imgaussfilt(noisy_scaled_fr,[.001 .5]) ;
    scaled_fr = avg_pref_fr .* max(real_psth(cc,:))./max(avg_pref_fr);
    cc_ind = sum(this_n(1:cc-1))+1:sum(this_n(1:cc));
    
    switch 1
        case 0
            sim_psth(cc,:) = mean(sim_rate_full(cc_ind,:));
            sim_rate = sim_rate_full;
        case 1
            [sim_raster, sim_rate, t_filt] = simulate_spikes_rates(this_fr, t, this_n(cc));
            sim_psth(cc,:) = mean(sim_rate);
            sim_rate_assembled(cc_ind,:) = sim_rate;
            sim_rate = sim_rate_assembled;
    end
    
end


s(1) = subplot(221)
plot(t, avg_pref_fr)
hold on
plot(t_filt,mean(sim_rate))
title('average firing rate')
ylabel('firing rate (spks/s)')
xlabel('time from cout (s)')
legend('real <FR>','simulated <FR> ','location','south')

s(2) = subplot(224)
[sim_psth_norm, ~, sim_pk_j] = sort_by_peak(sim_psth./max(sim_psth,[],2));
imagesc(sim_psth_norm,'x',t_filt)
title('simulated psth')
xlabel('time from cout (s)')
ylabel('cell # (sorted by peak time)')

s(3) = subplot(222)
[real_psth_norm, ~, real_pk_j] = sort_by_peak(real_psth./max(real_psth,[],2));
imagesc(real_psth_norm,'x',t)
title('real psth')
xlabel('time from cout (s)')
ylabel('cell # (sorted by peak time)')

%subplot(4,4,[9:10 13:14])
subplot(223)
plot(t_filt(sim_pk_j),t(real_pk_j),'.','markersize',20)
hold on
xlabel('simulated peak time (s)')
ylabel('real peak time (s)')
xlim(t([1 end]))
ylim(t([1 end]))
plot(xlim,xlim,'k')


% subplot(4,4,13:14)
% hist([t_filt(sim_pk_j); t(real_pk_j)]')


linkaxes(s,'x')

%colormap(flipud(bone))
colormap(hot)

suptitle(regions{rr})
%%
print(figure(12), fullfile('~/projects/ppcfof_tuning/sequences/',...
    [regions{rr} '_seq_v_simulated.jpg']),'-djpeg');
%%

%%


%%
%n_trials = round(n_trials./2);




figure(11); clf

plot(t,avg_pref_fr,'b','linewidth',2)
hold on
%plot(t,sim_fr,'r','linewidth',2)
plot(t_filt,nanmean(sim_raster_filt),'m','linewidth',2)

%%
sim_psth = nan(n_cells,length(t_filt));
for cc = 1:n_cells
    cc_ind = sum(n_trials(1:cc-1))+1:sum(n_trials(1:cc));
    sim_psth(cc,:) = mean(sim_raster_filt(cc_ind,:));
end

clf
%%
figure(12);
imagesc(sort_by_peak(sim_psth./max(sim_psth,[],2)),'x',t_filt)
colormap(hot)










%%
%ppc_cellids = cellids(region(:,1)=='p');

cells = cellids;
n_cells = length(cells);
taskp = zeros(length(cells),1);
lrp   = zeros(length(cells),1);
lrt   = zeros(length(cells),1);
rtbs = zeros(length(cells),1);
meanrate = zeros(length(cells),1);
%%

%%
ref_event = 'cpoke_end';
ycpk_all = cell(n_cells,1);
lcpk = cell(n_cells,1);
rcpk = cell(n_cells,1);
ystim = cell(n_cells,1);
lstim = cell(n_cells,1);
rstim = cell(n_cells,1);

peaktcpk = zeros(n_cells,1);
peaktstim = zeros(n_cells,1);

taskp_bl = zeros(n_cells,1);

fprintf('working on cell...')
tic
%%
for cc = 1256:1256%n_cells
    if mod(cc,5)==0
        fprintf('%i... ',cc)
        toc
        tic
    end
    tic
    [ycpk, tcpk, lcpk{cc}, rcpk{cc}, ~, ~] = make_pbups_raster(cells(cc),...
        'ref_event',ref_event,...
        'bin_size', bin_size, 'post', post, 'pre', pre);
    toc
    
    
    
    
    % separate psths for left/right, correct/error
    correct = vd.pokedR==(vd.bup_diff>0);
    left = vd.bup_diff<0;
    ylhit = ycpk(left&correct,:);
    yrhit = ycpk(~left&correct,:);
    ylerr = ycpk(~left&~correct,:);
    yrerr = ycpk(left&~correct,:);
    % find peak for each side correct trials
    [leftpk peakltcpk] = max(mean(ylhit));
    [rightpk peakrtcpk] = max(mean(yrhit));
    if leftpk>rightpk
        peaktcpk(cc) = peakltcpk;
    else
        peaktcpk(cc) = peakrtcpk;
    end
    
    
    
end


%% overall
figure(1); clf
plot(tcpk,mean(ycpk(vd.pokedR,:)))
hold on
plot(tcpk,mean(ycpk(~vd.pokedR,:)))

%%
sessn = bdata('select sessid from cells where cellid=7864');
sess_cells =  bdata('select cellid from cells where sessid={S}',sessn);

%sess_cells = cellids(1256:1257);
sess_cells =  cellids(region(:,1)=='f' & ...
    lr>.75 & mingamma<1.5); %& ~(normmean>1));

%%
unique_sessions = unique(sessids(region(:,1)=='f' & ...
    lr>.75 & mingamma<1.5));
%sess_cells = cellids(1:2);
fofperf = zeros(length(unique_sessions),1);
for ss = 1:length(unique_sessions)
    this_sess = unique_sessions(ss);
    %sd = get_sessdata(this_sess);
    peh = [sd.peh{1}.states];
    ncpks = cellfun(@rows, {peh.cpoke1});
    %ppccpkdur(ss) = diff(peh(find(ncpks==1,1)).cpoke1);
    fofperf(ss) = bdata('select total_correct from sessions where sessid={S}',this_sess);
end
%%
n_per_auc = .05/bin_size;

for ci = 1:2;%length(sess_cells)
    %%
    tic
    fprintf('working on cell %i of %i\n', ci, length(sess_cells));
    
    this_cellid = sess_cells(ci);
    [ycpk, tcpk, ~, ~, ~, vd] = make_pbups_raster(this_cellid,...
        'ref_event',ref_event,...
        'bin_size', bin_size, 'post', post, 'pre', pre);
    
    
    ypad =[ycpk zeros(size(ycpk,1), offset)];  % pad with extra zeros
    
    ycpk = filter(krn, 1, ypad, [], 2);
    ycpk = ycpk(:, 2*offset+1:end-1); % trim extra columns
    tcpk = tcpk(offset+1:end-1);
    ttotest = find(tcpk > -2 & tcpk < 1);
    auc_t = tcpk(1:n_per_auc:ttotest);
    
    
    if ci == 1
        auc    = zeros(length(sess_cells), length(ttotest)/n_per_auc);
        auc_p  = zeros(length(sess_cells), length(ttotest)/n_per_auc);
        auc_ci = zeros(length(sess_cells), length(ttotest)/n_per_auc);
        
        onseti = zeros(length(sess_cells), 1);
        lpsth   = zeros(length(sess_cells), length(ttotest));
        rpsth   = zeros(length(sess_cells), length(ttotest));
        prefpsth = zeros(length(sess_cells), length(ttotest));
        nonprefpsth = zeros(length(sess_cells), length(ttotest));
    end
    
    lpsth(ci, 1:length(ttotest)) = mean(ycpk(vd.pokedR==0,ttotest));
    rpsth(ci, 1:length(ttotest)) = mean(ycpk(vd.pokedR==1,ttotest));
    
    
    for ti = 1:length(ttotest)/n_per_auc
        ind = n_per_auc*(ti-1)+1:n_per_auc*(ti);
        tt = ttotest(ind);
        sc = sum(ycpk(:,tt),2);
        [this_auc, this_p, this_ci] = ...
            bootroc(sc(vd.pokedR==1),sc(vd.pokedR==0), 500, 95);
        auc(ci, ti) = this_auc;
        auc_p(ci, ti) = this_p;
        %auc_ci(ti) = this_ci;
        
    end
    
    
    if max(rpsth(ci,:)) > max(lpsth(ci,:))
        prefpsth(ci,1:length(ttotest)) = rpsth(ci,:);
        nonprefpsth(ci,1:length(ttotest)) = lpsth(ci,:);
    else
        prefpsth(ci,1:length(ttotest)) = lpsth(ci,:);
        nonprefpsth(ci,1:length(ttotest)) = rpsth(ci,:);
    end
    toc
end
%%
all_cellids = cellids;
%%


%% significant auc plot
align_ind = 8;
[align_strs, align_args] = dyn_align_LUT;
align_str = align_strs{align_ind};

f = load(fullfile('~/projects/ppcfof_tuning/sequences/',...
    ['all_fof_aucs_psths_align_' align_str '.mat']));

p = load(fullfile('~/projects/ppcfof_tuning/sequences/',...
    ['all_ppc_aucs_psths_align_' align_str '.mat']));

s = load(fullfile('~/projects/ppcfof_tuning/sequences/',...
    ['all_striatum_aucs_psths_align_' align_str '.mat']));


minmeanfr = 1;
f.ind = mean(f.pref_psth,2) > minmeanfr;
p.ind = mean(p.pref_psth,2) > minmeanfr;
s.ind = mean(s.pref_psth,2) > minmeanfr;

[f.phat  f.pci] = binofit(sum(abs(f.auc_p(f.ind,:)-.5)>.499,1),...
    repmat(size(f.auc(f.ind,:),1),1,size(f.auc(f.ind,:),2)))
[p.phat  p.pci] = binofit(sum(abs(p.auc_p(p.ind,:)-.5)>.499,1),...
    repmat(size(p.auc(p.ind,:),1),1,size(p.auc(p.ind,:),2)))
[s.phat  s.pci] = binofit(sum(abs(s.auc_p(s.ind,:)-.5)>.499,1),...
    repmat(size(s.auc(s.ind,:),1),1,size(s.auc(s.ind,:),2)))

fh = figure(102); clf
plot(f.t,f.phat,'r','linewidth',2)

hold on
plot(p.t,p.phat,'k','linewidth',2)
plot(s.t,s.phat,'b','linewidth',2)
plot(p.t,p.pci(:,1),':k','linewidth',2)
plot(p.t,p.pci(:,2),':k','linewidth',2)
plot(f.t,f.pci(:,1),':r','linewidth',2)
plot(f.t,f.pci(:,2),':r','linewidth',2)
plot(s.t,s.pci(:,1),':b','linewidth',2)
plot(s.t,s.pci(:,2),':b','linewidth',2)

legend('fof','ppc','striatum')
ylabel('% cells significantly encoding choice')
xlabel(align_str)
print(fh,['percent_sign_cells_' align_str '.pdf'],'-dpdf','-bestfit')
%%

sig_cells = cellids;

d= p;

titlestr = 'PPC significant peaks';
savename = 'ppc_sigpeaks_auc_vs_psth'
d.sig_cells = ismember(d.cells,sig_cells);

which_cells = (d.sig_cells);

figure(1); clf

prefpsth = d.prefpsth(which_cells,:);
nonprefpsth = d.nonprefpsth(which_cells,:);
lpsth = d.lpsth(which_cells,:);
rpsth = d.rpsth(which_cells,:);
auc_p    = d.auc_p(which_cells, :);
auc      = d.auc(which_cells, :);
tcpk     = d.tcpk;
ttotest  = d.ttotest;
bin_size = diff(d.tcpk(1:2));


toplot = mean(prefpsth(:,tcpk>-1.5 & tcpk<0),2)>2;
cumsig = cumsum(abs(auc_p-.5)>.499,2)==ceil(.15/bin_size);
cumsig(:,end)=ceil(.15/bin_size);
[i, j] = find(cumsum(cumsum(cumsig,2),2)==1);

colormap(hot)

i = i(ismember(i,find(toplot)));


subplot(231)
imagesc(abs(auc(i,:)-.5),'x',tcpk(ttotest),[0 .5]);
hold on;plot([0 0], ylim, 'w','linewidth',2)
ylabel('cell # (sorted by selectivity onset)')
title('side-selectivity')
cb=colorbar
title(cb, '|auc|-.5')

subplot(232)
imagesc(prefpsth(i,:)./max(prefpsth(i,:),[],2),'x',tcpk(ttotest),[0 1]);
hold on;plot([0 0], ylim, 'w','linewidth',2)
title('preferred side psth')
cb=colorbar
title(cb, 'norm fr')

subplot(233)
imagesc(nonprefpsth(i,:)./max(prefpsth(i,:),[],2),'x',tcpk(ttotest),[0 1]); colorbar
hold on;plot([0 0], ylim, 'w','linewidth',2)
title('non-preferred side psth')
cb=colorbar
title(cb, 'norm fr')

[~, maxi] = max(prefpsth,[],2);
[~, i]    = sort(maxi);
subs = sub2ind(size(prefpsth), 1:length(maxi), maxi');
subs0 = sub2ind(size(prefpsth), 1:length(maxi), max(1,-2+maxi'));
subs1 = sub2ind(size(prefpsth), 1:length(maxi), max(1,-1+maxi'));
subs2 = sub2ind(size(prefpsth), 1:length(maxi), min(size(prefpsth,2),1+maxi'));
subs3 = sub2ind(size(prefpsth), 1:length(maxi), min(size(prefpsth,2),2+maxi'));


toplot = (abs(auc_p(subs)-.5) > .499) + ...
    (abs(auc_p(subs0)-.5) > .499) + ...
    (abs(auc_p(subs1)-.5) > .499) + ...
    (abs(auc_p(subs2)-.5) > .499) + ...
    (abs(auc_p(subs3)-.5) > .499)== 5;% & ;



%toplot = toplot & d.tcpk(d.ttotest(maxi)) < 0;



toplot = toplot & mean(prefpsth(:,tcpk(d.ttotest)>-2 & tcpk(d.ttotest)<2),2)'>1;

i = i(ismember(i,find(toplot)));
%plot(maxi(i),1:length(maxi),'.')

subplot(234)
imagesc(abs(auc(i,:)-.5),'x',tcpk(ttotest),[0 .5]);
hold on;plot([0 0], ylim, 'w','linewidth',2)
xlabel('time from cout (s)')
ylabel('cell # (sorted by preferred side peak)')
cb=colorbar
title(cb, '|auc|-.5')

subplot(235)
imagesc(prefpsth(i,:)./max(prefpsth(i,:),[],2),'x',tcpk(ttotest),[0 1]);
hold on;plot([0 0], ylim, 'w','linewidth',2)
xlabel('time from cout (s)')
cb=colorbar
title(cb, 'norm fr')

subplot(236)
imagesc(nonprefpsth(i,:)./max(prefpsth(i,:),[],2),'x',tcpk(ttotest),[0 1]);
hold on;plot([0 0], ylim, 'w','linewidth',2)
xlabel('time from cout (s)')
cb=colorbar
title(cb, 'norm fr')

colormap(hot.^.5)

suptitle(titlestr)
print(figure(1), ['~/projects/ppcfof_tuning/sequences/' savename '.jpg'],'-djpeg')
%%

figure(1); clf

%%
sig_cells = 7752;

d= f;
subplot(121)
cla

titlestr = 'FOF';
savename = 'fof_sigpeaks_stimper_auc_vs_psth'
d.sig_cells = ismember(d.cells,sig_cells);

which_cells = (d.sig_cells);

prefpsth = d.prefpsth(which_cells,:);
nonprefpsth = d.nonprefpsth(which_cells,:);
lpsth = d.lpsth(which_cells,:);
rpsth = d.rpsth(which_cells,:);
auc_p    = d.auc_p(which_cells, :);
auc      = d.auc(which_cells, :);
tcpk     = d.tcpk;
ttotest  = d.ttotest;
bin_size = diff(d.tcpk(1:2));


[~, maxi] = max(prefpsth,[],2);
[~, i]    = sort(maxi);
subs = sub2ind(size(prefpsth), 1:length(maxi), maxi');
subs0 = sub2ind(size(prefpsth), 1:length(maxi), max(1,-2+maxi'));
subs1 = sub2ind(size(prefpsth), 1:length(maxi), max(1,-1+maxi'));
subs2 = sub2ind(size(prefpsth), 1:length(maxi), min(size(prefpsth,2),1+maxi'));
subs3 = sub2ind(size(prefpsth), 1:length(maxi), min(size(prefpsth,2),2+maxi'));


toplot = (abs(auc_p(subs)-.5) > .499) + ...
    (abs(auc_p(subs0)-.5) > .499) + ...
    (abs(auc_p(subs1)-.5) > .499) + ...
    (abs(auc_p(subs2)-.5) > .499) + ...
    (abs(auc_p(subs3)-.5) > .499)== 5;% & ;


toplot = toplot & d.tcpk(d.ttotest(maxi)) < 0;

toplot = toplot & mean(prefpsth(:,tcpk>-1.5 & tcpk<0),2)'>1;


i = i(ismember(i,find(toplot)));
%plot(maxi(i),1:length(maxi),'.')



plot(d.tcpk(d.ttotest),(length(i):-1:1)+2*prefpsth(i,:)'./max(prefpsth(i,:),[],2)','linewidth',2,'color','g')
hold on
plot(d.tcpk(d.ttotest),(length(i):-1:1)+2*nonprefpsth(i,:)'./max(prefpsth(i,:),[],2)','linewidth',2,'color','r')
plot([0 0],ylim,'k')
colormap(hot.^.5)
%xlim([-3 3])
title(titlestr)
xlabel('time from cout (s)')
ylabel('cell # (sorted by preferred peak')
%print(figure(1), ['~/projects/ppcfof_tuning/sequences/' 'selective_sequence_during_stim_linepsth' '.jpg'],'-djpeg')




%%
tic
[ystim{cc}, tstim, lstim{cc}, rstim{cc}, ~, ~] = make_pbups_raster(cells(cc),...
    'ref_event',ref_event,...
    'bin_size', bin_size, 'post', post, 'pre', pre);
toc
ycpk_all{cc} = ycpk;


[rtb, p, ~] = test_ridge_to_bg(ycpk,1000);
taskp_bl(cc) = p;

%%

fprintf('working on cell...')
for cc = 1:length(cells)
    fprintf('%i... ',cc)
    [ycpk, tcpk, lcpk, rcpk, ad, vd] = make_pbups_raster(cells(cc),'ref_event','cpoke_end',...
        'bin_size', bin_size, 'post', post, 'pre', pre);
    correct = vd.pokedR==(vd.bup_diff>0);
    left = vd.bup_diff<0;
    
    ylhit = ycpk(left&correct,:);
    [ylhittdur, ylhitsortorder] = sort(vd.T(left&correct));
    
    yrhit = ycpk(~left&correct,:);
    [yrhittdur, yrhitsortorder] = sort(vd.T(~left&correct));
    
    ylerr = ycpk(~left&~correct,:);
    [ylerrtdur, ylerrsortorder] = sort(vd.T(~left&~correct));
    
    yrerr = ycpk(left&~correct,:);
    [yrerrtdur, yrerrsortorder] = sort(vd.T(left&~correct));
    
    leftpk = max(mean(ylhit));
    rightpk = max(mean(yrhit));
    if leftpk>rightpk
        [rtb, p, ~] = test_ridge_to_bg(ylhit,1000);
    else
        [rtb, p, ~] = test_ridge_to_bg(yrhit,1000);
    end
    taskp(cc) = p;
    rtbs(cc)  = rtb;
    
    tind = tf<0&tf>-.5;
    
    leftct = sum(ylhit(:,tind),2);
    rightct = sum(yrhit(:,tind),2);
    [~, lrp(cc), ~, stats] = ttest2(leftct,rightct);
    lrt(cc) = stats.tstat;
    
    meanrate(cc) = mean(ycpk(:));
    
    evidence = (cumsum(rcpk,2)-cumsum(lcpk,2));
    
    
    ypad =[ycpk zeros(size(ycpk,1), offset)];  % pad with extra zeros
    
    yf = filter(krn, 1, ypad, [], 2);
    yf = yf(:, 2*offset+1:end-1); % trim extra columns
    tf = tcpk(offset+1:end-1);
    evf = evidence(:,offset+1:end-1);
    
    %%
    f = figure(1); clf
    subplot(321)
    
    evplot = evf(:,tind);
    yplot = yf(:,tind);
    tfplot = repmat(tf(tind),size(yplot,1),1) ;
    trialfplot = repmat((1:size(yplot,1))',1,size(yplot,2)) ;
    scatter(evplot(:)-.125+rand(size(evplot(:)))/4, yplot(:),[],trialfplot(:))
    colormap(flipud(colormapRedBlue.*[.8 .8 .8]))
    colorbar
    xlim(max(abs(evplot(:))).*[-1 1])
    title('evidence tuning')
    ylabel('spk/s')
    xlabel('click difference (R-L)')
    
    evind = evplot(:)-min(evplot(:))+1;
    
    subplot(322)
    m=accumarray(evind, yplot(:), [], @mean)
    msem=accumarray(evind, yplot(:), [], @sem)
    pev = accumarray(evind, ones(size(yplot(:))), [], @sum)
    pev = pev ./ sum(pev)
    
    mx = min(evplot(:))+[0:length(m)-1];
    plot(mx, m, '-')
    hold on
    plot(mx, m+msem,'-k')
    plot(mx, m-msem,'-k')
    title('evidence tuning')
    ylabel('spks/s')
    xlabel('click difference (R-L)')
    xlim(max(abs(evplot(:))).*[-1 1])
    % correct trials
    
    subplot(323)
    imagesc([ylhit(ylhitsortorder,:); -yrhit(yrhitsortorder,:)],'x',tcpk)
    hold on
    scatter(-ylhittdur, 1:length(ylhittdur), '.k')
    scatter(-yrhittdur, length(ylhittdur)+(1:length(yrhittdur)), '.k')
    caxis([-2 2])
    
    %colormap(flipud(bone))
    hold on
    plot(xlim, [1 1].*size(ylhit,1),'k')
    title('raster (error trials)')
    ylabel('trials sorted by choice')
    xlabel('time from cpoke end (s)')
    
    
    subplot(325)
    
    mylf = nanmean(yf(left&correct,:));
    semylf = nansem(yf(left&correct,:));
    myrf = nanmean(yf(~left&correct,:));
    semyrf = nansem(yf(~left&correct,:));
    
    plot(tf, mylf,'b')
    hold on
    plot(tf, myrf,'r')
    plot(tf, mylf+semylf,'b')
    plot(tf, mylf-semylf,'b')
    
    
    
    plot(tf, myrf+semyrf,'r')
    plot(tf, myrf-semyrf,'r')
    
    title('PSTH (correct trials)')
    xlabel('time from cpoke end (s)')
    ylabel('spks/s')
    legend({'went left','went right'})
    plot([0 0], ylim, 'k')
    % plot inccorrect trials
    
    
    subplot(324)
    imagesc([ylerr; -yrerr],'x',tcpk)
    hold on
    scatter(-ylerrtdur, 1:length(ylerrtdur), '.k')
    scatter(-yrerrtdur, length(ylerrtdur)+(1:length(yrerrtdur)), '.k')
    caxis([-2 2])
    colormap(flipud(colormapRedBlue))
    
    hold on
    plot(xlim, [1 1].*size(ylerr,1),'k')
    title('raster (error trials)')
    ylabel('trials sorted by choice')
    xlabel('time from cpoke end (s)')
    
    subplot(326)
    mylf = nanmean(yf(~left&~correct,:));
    semylf = nansem(yf(~left&~correct,:));
    myrf = nanmean(yf(left&~correct,:));
    semyrf = nansem(yf(left&~correct,:));
    
    plot(tf, mylf,'--b')
    hold on
    
    plot(tf, myrf,'--r')
    plot(tf, mylf+semylf,'--b')
    plot(tf, mylf-semylf,'--b')
    
    
    
    plot(tf, myrf+semyrf,'--r')
    plot(tf, myrf-semyrf,'--r')
    pause(.1)
    
    plot([0 0], ylim, 'k')
    title('PSTH (error trials)')
    xlabel('time from cpoke end (s)')
    ylabel('spks/s')
    legend({'went left','went right'})
    
    linkaxes([subplot(323) subplot(324) subplot(325) subplot(326)],'x')
    
    assert(length(region)==length(cells))
    
    if region(cc,1) == 'p'
        st = ['PPC ' num2str(cells(cc))];
    elseif region(cc,1) == 'f'
        st = ['FOF ' num2str(cells(cc))];
    elseif region(cc,1) == 'd'
        st = ['dH ' num2str(cells(cc))];
    end
    
    suptitle(st);
    
    print(f, ['/Users/oroville/projects/ppcfof_tuning/figures/evidencepsths/'...
        strrep(st,' ','_') '.pdf'] , '-dpdf', '-bestfit')
    %%
    tbin = repmat(1:sum(tind),size(yplot,1),1);
    mt = accumarray(tbin(:), yplot(:), [], @mean);
    pt = accumarray(tbin(:), ones(size(yplot(:))), [], @sum);
    pt = pt./sum(pt);
    %skaggs evidence information
    evinfo = nansum(m .* log2(m ./ mean(yplot(:))) .* pev)
    %skaggs evidence information
    tinfo = nansum(mt .* log2(mt ./ mean(yplot(:))) .* pt)
    
    
end











%% load psths for sequence plot

% hmmm something weird happens if you use cpokeout rather than
% cpokeout-nomask

figure(1); clf
ax =axes();

opts.alignment   = 'cpokeout-nomask';
%opts.alignment   = 'stimstart-cout-mask';
opts.flip_bool   = 'false';
opts.select_bool = 'data.trials.hit==1';%& data.trials.correct_dir~=0';
opts.norm_type = 'none';

opts.repack = false;


opts.condsort    = 'data.trials.correct_dir';
opts.bin_conds = false;

[~, res_cout] = pbups_psth(cellids,'condsort', opts.condsort,...
    'select_bool', opts.select_bool,...
    'flip_bool',opts.flip_bool,'norm_type',opts.norm_type,...
    'bin_conds',opts.bin_conds,'axes_in',ax, 'alignment', opts.alignment, 'repack', opts.repack);
%%
opts.alignment   = 'stimstart2';
[~, res_stim] = pbups_psth(cellids,'condsort', opts.condsort,...
    'select_bool', opts.select_bool,...
    'flip_bool',opts.flip_bool,'norm_type',opts.norm_type,...
    'bin_conds',opts.bin_conds,'axes_in',ax, 'alignment', opts.alignment, 'repack', opts.repack);
%%
opts.alignment   = 'cpokein';
[~, res_cin] = pbups_psth(cellids,'condsort', opts.condsort,...
    'select_bool', opts.select_bool,...
    'flip_bool',opts.flip_bool,'norm_type',opts.norm_type,...
    'bin_conds',opts.bin_conds,'axes_in',ax, 'alignment', opts.alignment, 'repack', opts.repack);
%%
opts.alignment   = 'spokein';
[~, res_sin] = pbups_psth(cellids,'condsort', opts.condsort,...
    'select_bool', opts.select_bool,...
    'flip_bool',opts.flip_bool,'norm_type',opts.norm_type,...
    'bin_conds',opts.bin_conds,'axes_in',ax, 'alignment', opts.alignment, 'repack', opts.repack);

%% make sequence plot
%dirs = {'left', 'right'};

align_sort = 'cout';

align_sort = 'cin';
align_sort = 'sin';


align_plot = align_sort;

splitsort  = true;
group_pref = true;
sig_only   = true;

ppc = region(:,1)=='p';
fof = region(:,1)=='f';

region_name = 'fof';

res_sort = eval(['res_' align_sort]);
res_plot = eval(['res_' align_plot]);


savename = ['lr_' region_name '_' align_plot 'sortedby_' align_sort];


tind_plot = res_plot.x>-4 & res_plot.x<3;
tind_sort = res_sort.x>-10;


h = figure(1); clf

region_ind = eval(region_name);
plot_ind = region_ind & taskmod ;%& mingamma<1.5;
if sig_only
    savename = [ 'sig' savename];
    plot_ind = plot_ind & prefp<.05 & taskmod ;%& mingamma<1.5;
end


set(0,'defaultaxesfontsize',15)


cmap = flipud(colormapBlues).^.65;

if ~splitsort
    sortfun = @(psth) squeeze(max(psth(plot_ind,:,tind_sort),[],2));
    
else
    %savename = ['lr_' align_plot 'sortedby_side' align_sort];
    sortfun = @(psth) permute(psth(plot_ind,:,tind_sort),[1 3 2]);
end

sorting_psth = sortfun(res_sort.y_mean_eachCell);
plot_psth = permute(res_plot.y_mean_eachCell(plot_ind,:,tind_plot),[1 3 2]);

ppc_psth_flat = reshape(sorting_psth,[size(sorting_psth,1) ...
    size(sorting_psth,2)*size(sorting_psth,3)]);
[pk, pki]  = max(horzcat(ppc_psth_flat),[],2);

if group_pref
    m = size(sorting_psth,2);
    prefside = pki <= m;
    plot_psth(~prefside,:,:) = plot_psth(~prefside,:, [2 1]);
    [pk, pki]  = max(plot_psth(:,:,1),[],2);
    cond_labels = {'preferred correct trials','non-preferred correct trials'}
end

trep = repmat(res_sort.x,size(plot_psth,1),1);

t_comind = (sum(trep.*plot_psth(:,:,1)./sum(plot_psth(:,:,1),2),2))

%t_com = res_sort.x(t_comind);



%[pk, pki]  = max(horzcat(plot_psth),[],2);
sortbyt = 'peak';
switch sortbyt
    case 'peak'
        
        [~, sortix] = sort(pki);
    case 'com'
        [~, sortix] = sort(t_comind);
end

n_conds = size(plot_psth,3);

for cc = 1:n_conds
    s(cc) = subplot(1,n_conds,cc)
    imagesc(plot_psth(sortix,:,cc)./pk(sortix,:),'x',res_plot.x(tind_plot))
    
    xlabel(['time from ' align_plot ' (s)'])
    if cc==1
        ylabel(['cell # (sorted by peak firing time from ' align_sort ')'])
    end
    title(cond_labels{cc},'fontsize',12)
    cb = colorbar
    title(cb,'norm fr')
    caxis([0 1])
    
    hold on
    plot([0 0],ylim,'k','linewidth',2)
    
end



linkaxes(s)

colormap(cmap)
%suptitle([upper(region_name)])

%this_savename = [savename];

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

print(figure(1), fullfile('~/projects/ppcfof_tuning',savename),'-dpdf')






figure(2); clf


xx = [res_plot.x res_plot.x]

histogram(xx(pki),25,'normalization','probability')
title([region_name ' peak times'])
xlabel('time from cpoke end (s)')
ylabel('fraction of cells')

histsavename = [this_savename '_peakhist'];


print(figure(2), fullfile('~/projects/ppcfof_tuning',histsavename),'-dpdf')







%%
%opts.alignment = 'stimstart'
%opts.alignment = 'cpokeend'
opts.alignment = 'cpokeout-nomask'
opts.bin_conds = true;
opts.cond_edges = [-28 -21 -14 -7 0 7 14 21 28 ];
opts.select_bool = 'data.trials.hit==1';
[~, resh] = pbups_psth(cellids,'condsort', opts.condsort,...
    'select_bool', opts.select_bool,...
    'flip_bool',opts.flip_bool,'norm_type',opts.norm_type,...
    'bin_conds',opts.bin_conds,'cond_edges',opts.cond_edges,...
    'condsort','data.trials.bup_diff',...
    'axes_in',ax, 'alignment', opts.alignment, 'repack', opts.repack);


opts.select_bool = 'data.trials.hit==0';

[~, rese] = pbups_psth(cellids,'condsort', opts.condsort,...
    'select_bool', opts.select_bool,...
    'flip_bool',opts.flip_bool,'norm_type',opts.norm_type,...
    'bin_conds',opts.bin_conds,'cond_edges',opts.cond_edges,...
    'condsort','data.trials.bup_diff',...
    'axes_in',ax, 'alignment', opts.alignment, 'repack', opts.repack);


cmin = 0;
cmax = .6;
pwr = .6;
n_conds =  size(resh.y_mean_eachCell,2);
grad = linspace(cmin,cmax,n_conds/2)'.^pwr;
gradF = flipud(grad);
clrs = [grad grad ones(size(grad));  ones(size(grad)) gradF gradF];
% should save this thing in ppcfoftuning
% savename = [cel'lr_' opts.alignment '_psth'];
%
% save(fullfile('~/projects/ppcfof_tuning',[savename '.mat']),'opts','res')


plot_psth = permute(resh.y_mean_eachCell(:,:,:),[1 3 2]);
% plot psths by gamma for all cells




% grad2 = linspace(cmin,cmax,size(clrs2,1)/2)';
% gradF2 = flipud(grad2);
% clrs2 = [grad2 grad2 ones(size(grad2));  ones(size(grad2)) gradF2 gradF2];

figure(2); clf


%sig_select_str = ['lr>0.75 & mingamma<=1.5 & prefp_stimcout<0.05 & normmean>1'];
%sig_cells = cell2mat(extracting(cell_list, 'cellid', sig_select_str));


% for cc = 1:length(these_cellids)
%     cind = find(cellids==these_cellids(cc));
for cc = 1;%1:length(cellids)
    %     is_sig = ismember(cellids(cc), sig_cells)
    %     if is_sig
    %         sigstr = 'sig';
    %     else
    %         sigstr= '';
    %     end
    %
    %     cind = cc;
    for dd = 1:size(plot_psth,3)
        plot(resh.x,plot_psth(cind,:,dd),'color',clrs(dd,:),'linewidth',2)
        hold on
    end
    % title(sprintf('%s %i, prefp stimcout=%.3f %s' ,region(cc,:),cellids(cc),prefp(cc),sigstr ))
    plot([ 0 0], ylim,'k')
    plot([ -1.2 -1.2], ylim,'k')
    this_savename = sprintf('%s%s%i_bygamma_%s.jpg' ,'sig','ppc',1755,opts.alignment )
    %this_savename = sprintf('%s%s%i_bygamma.jpg' ,sigstr,region(cc,:),cellids(cc) )
    ylabel('spks/s')
    xlabel('time from center off (s)')
    print(figure(2), fullfile('~/projects/ppcfof_tuning',this_savename),'-djpeg')
    
    
end

%% plot "tuning curve" by gamma
hit_psth = permute(resh.y_mean_eachCell(:,:,:),[1 3 2]);
err_psth = permute(rese.y_mean_eachCell(:,:,:),[1 3 2]);


cells_toplot = find((hitprefp < .05 | prefp<0.05) & mingamma < 1.5);


ctrs = mean([opts.cond_edges(2:end); opts.cond_edges(1:end-1)]);
ctrs_flip = fliplr(ctrs);
clrs_flip = flipud(clrs);
tind = res.x>-.3 & res.x<.1;
for cc = 1:length(cells_toplot)
    cind = cells_toplot(cc);
    is_sig = ismember(cellids(cind), sig_cells)
    if is_sig
        sigstr = 'sig';
    else
        sigstr= '';
    end
    figure(3); clf
    
    plot(ctrs,squeeze(mean(hit_psth(cind,tind,:))),...
        '-','linewidth',2, 'color',[1 1 1].*.75)
    hold on
    plot(ctrs,squeeze(mean(err_psth(cind,tind,:))),...
        '--','linewidth',2, 'color',[1 1 1].*.75)
    h = legend('hits','errors','location','eastoutside')
    
    for dd = 1:size(hit_psth,3)
        plot(ctrs(dd),mean(hit_psth(cind,tind,dd)),'.','markersize',30,'color',clrs(dd,:))
        plot(ctrs(dd),mean(err_psth(cind,tind,dd)),'.','markersize',30,'color',clrs(dd,:))
    end
    title(sprintf('%s %i' ,region(cind,:),cellids(cind) ))
    xlabel('gamma bin')
    ylabel('firing rate')
    this_savename = sprintf('%s%i_meanbygamma.jpg' ,region(cind,:),cellids(cind) )
    print(figure(3), fullfile('~/projects/ppcfof_tuning',this_savename),'-djpeg')
    pause(1)
end












%%
opts.select_bool = 'data.trials.hit==1 & data.trials.correct_dir~=0';
[~, resh] = pbups_psth(cellids,'condsort', opts.condsort,...
    'select_bool', opts.select_bool,...
    'flip_bool',opts.flip_bool,'norm_type',opts.norm_type,...
    'bin_conds',opts.bin_conds,'axes_in',ax, 'alignment', opts.alignment);
opts.select_bool = 'data.trials.hit==0 & data.trials.correct_dir~=0';
[~, rese] = pbups_psth(cellids,'condsort', opts.condsort,...
    'select_bool', opts.select_bool,...
    'flip_bool',opts.flip_bool,'norm_type',opts.norm_type,...
    'bin_conds',opts.bin_conds,'axes_in',ax, 'alignment', opts.alignment);


%%
if 0
    switch 0
        case 0
            savename = 'lr_cpokeout_psth';
            cond_labels = {'correct left trials','correct right trials'}
        case 1
            savename = 'lr_easyhard_cpokeout_psth';
            cond_labels = {'easy left trials','hard left trials',...
                'hard right trials','easy right trials'};
    end
    
    load(fullfile('~/projects/ppcfof_tuning',[savename '.mat']),'opts','res')
end
%%
dirs = {'left', 'right'};
%cond_labels = {'eror left trials','error right trials'}

h = figure(1); clf

ppc = region(:,1)=='p';
fof = region(:,1)=='f';

region_name = 'ppc';
plot_ind = eval(region_name);

set(0,'defaultaxesfontsize',15)


cmap = flipud(colormapBlues).^.65;


psthh = permute(resh.y_mean_eachCell(plot_ind,:,:),[1 3 2]);
psthe = permute(rese.y_mean_eachCell(plot_ind,:,:),[1 3 2]);
sorting_psth = [psthh ];
ppc_psth_flat = reshape(sorting_psth,[size(sorting_psth,1) size(sorting_psth,2)*size(sorting_psth,3)])

[pk, pki]  = max(horzcat(ppc_psth_flat),[],2);

[~, sortix] = sort(pki);

n_conds = size(sorting_psth,3);


for cc = 1:n_conds
    subplot(2,n_conds,cc)
    imagesc(psthh(sortix,:,cc)./pk(sortix,:),'x',res.x)
    xlabel('time from center out (s)')
    if cc==1
        ylabel('cell # (sorted by peak firing time)')
    end
    %title(cond_labels{cc},'fontsize',12)
    title(['correct ' dirs{cc} ' trials']);
    cb = colorbar
    title(cb,'norm fr')
    caxis([0 1])
end

for cc = 1:n_conds
    subplot(2,n_conds,2+cc)
    this_cond = n_conds-cc+1;
    imagesc(psthe(sortix,:,this_cond)./pk(sortix,:),'x',res.x)
    xlabel('time from center out (s)')
    if cc==1
        ylabel('cell # (sorted by peak firing time)')
    end
    %title(cond_labels{cc},'fontsize',12)
    title(['error ' dirs{cc} ' trials']);
    cb = colorbar
    title(cb,'norm fr')
    caxis([0 1])
end

colormap(cmap)
%suptitle([upper(region_name)])

this_savename = [region_name '_' savename];

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

print(figure(1), fullfile('~/projects/ppcfof_tuning',this_savename),'-dpdf')

%%
subplot(141)
imagesc(psth1(sortix,:)./pk(sortix,:),'x',res.x)
colorbar
title('PPC')

subplot(142)
imagesc(psth2(sortix,:)./pk(sortix,:),'x',res.x)
colorbar
colormap(cmap)

title('PPC')

%figure(2); clf

psth1 = squeeze(res.y_mean_eachCell(fof,1,:));
psth2 = squeeze(res.y_mean_eachCell(fof,2,:));

[pk, pki]  = max(horzcat(psth1,psth2),[],2);

[~, sortix] = sort(pki);

subplot(143)
imagesc(psth1(sortix,:)./pk(sortix,:),'x',res.x)
colorbar
title('FOF')

subplot(144)
imagesc(psth2(sortix,:)./pk(sortix,:),'x',res.x)
colorbar
colormap(cmap)

title('FOF')

