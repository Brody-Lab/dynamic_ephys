close all
clear all
dp = set_dyn_path

%% get full cell_list; *all* single units recorded in pbups project
cell_list = dyn_cells_db;   % That has to be run once to create cell_list

%% compile a bunch of psths aligned to the center poke time and the center out time
% produce a stim on aligned plot and a center out aligned plot
align_strs      = dyn_align_LUT;
cout_align_ind  = find(ismember(align_strs,'cpokeout'));
cin_align_ind   = find(ismember(align_strs,'stimstart-cout-mask'));


% turn long_trials_only on if you want to only use the long trials for the
% figures
long_trial_dur = 1;
long_trials_only = false;

select_str = 'strcmp(region,''fof'') & normmean > .5' ;
%select_str = 'cellid==18181'
cellids = cell2mat(extracting(cell_list, 'cellid', select_str));
psr = cell2mat(extracting(cell_list, 'prefsideright', select_str));
prefp = cell2mat(extracting(cell_list, 'prefp', select_str));
normmean = cell2mat(extracting(cell_list, 'normmean', select_str));
ncells = size(cellids);
%%
rind = 3;
lind = 2;

for cc = 1:ncells
    try
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

%%
%%%% PLOT SEQUENCE PLOT
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

%% PREF/NONPREF pop average PSTH
mn_fr       = nanmean(cin_psth_hit(:,:,1,1),2);
good_cint   = cin_t >= -.5 & cin_t <= 1.5;
good_coutt  = cout_t >= -1.5 & cout_t <= 1;

switch 1
    case 0
        good_cells  = mn_fr > 1 & prefp < .01;
        pref_r = nanmean(cin_psth_hit(:,:,rind,1),2) > nanmean(cin_psth_hit(:,:,lind,1),2);
        pref_l = ~pref_r;
    case 1
        good_cells  = mn_fr > 1 & prefp < .01;
        pref_r = psr == 1;
        pref_l = 0   == psr;
           
end

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


pref_color =[.5 .3 .75];
npref_color = hsv2rgb(rgb2hsv(1-pref_color).*[1 .25 1]);


fh = figure(2); clf

ax(1) = subplot(121);hold(ax(1),'on');
ax(2) = subplot(122);hold(ax(2),'on');

plot(ax(1),[ 0 0], [0 100],'k')
plot(ax(2),[ 0 0], [0 100],'k')
ylims = ([floor(min(psths(:))*10)/10 ceil(max(psths(:))*10)/10])
ylim(ax(1),ylims)
ylim(ax(2),ylims)
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



cin_pref_psth_good = cin_pref_psth;
cin_npref_psth_good = cin_npref_psth;

cout_pref_psth_good = cout_pref_psth;
cout_npref_psth_good = cout_npref_psth;

pref_combo = [cin_pref_psth_good cout_pref_psth_good];

normsort = @(A,B) sort_by_peak(norm_by_peak(A,B),B);
%normsort = @(A,B) sort_by_peak(A,B);
cax = [0 1]

fh=figure(3); clf
set(fh,'position',[2 5 4.5 6 ],'papersize',[5 3],'paperpositionmode','auto')

s(1)=subplot(221)
imagesc(normsort(cin_pref_psth_good, pref_combo),'x',cin_t(good_cint))
caxis(cax)
s(2)=subplot(222)
imagesc(normsort(cout_pref_psth_good, pref_combo),'x',cout_t(good_coutt))
caxis(cax)

s(3)=subplot(223)
imagesc(normsort(cin_npref_psth_good, pref_combo),'x',cin_t(good_cint))
caxis(cax)
s(4)=subplot(224)
imagesc(normsort(cout_npref_psth_good, pref_combo),'x',cout_t(good_coutt))
caxis(cax)
colormap(flipud(gray.^.7))
ylabel(subplot(221),'cell # (sorted by peak)')
ylabel(subplot(223),'cell # (sorted as above)')
hold(s(1),'on')
hold(s(2),'on')
hold(s(3),'on')
hold(s(4),'on')
plot(s(1),[ 0 0], [0 1000],'k')
plot(s(2),[ 0 0], [0 1000],'k')
plot(s(3),[ 0 0], [0 1000],'k')
plot(s(4),[ 0 0], [0 1000],'k')

set(subplot(221),'ytick',[300 600])
set(subplot(223),'ytick',[300 600])
set(subplot(222),'ytick',[])
set(subplot(224),'ytick',[])
set(subplot(222),'ytick',[])
set(subplot(224),'ytick',[]) 
xlabel(subplot(223),'time from stim onset (s)')
xlabel(subplot(224),'from movement (s)')

cb = colorbar
cb.Position = cb.Position + [0 .3 .052 -.2]

%colormap(flipud(colormapBlues.^.5))
%%
good_cells  = normmean > 5 & prefp < .01;

edges = [-2 -.8:.2:.8 2];
pref_color  = [.8 .25 .8];
npref_color = [.8 .65 .25];
example_cell_psth('cells',cellids(good_cells),'meta',1,'type','chrono',...
    'edges',edges,'norm','onset','top_color',pref_color,'bot_color',npref_color)

%% try dprime plot

dp_cin = nan(size(cin_psth(:,:,1,1)));
dp_cout = nan(size(cout_psth(:,:,1,1)));

cc = find(prefp<.0001 & good_cells,1);

good_cell_ix = find(good_cells);
good_cell_ix = 1:ncells;
ngood = length(good_cell_ix);

nboot = 250;

cout_auc = nan(ngood,nt);
cout_p   = nan(ngood,nt);
cout_ci  = nan(ngood,nt,2);
fprintf('working on cell...')
tic
for cc = 1:ngood
    if mod(cc,5)==0
        toc;
        fprintf('%i...',cc);
        tic;
    end
    this_id = cellids(good_cell_ix(cc));
    
    d           = dyn_cell_packager(this_id);
    cin_frates  = d.frate{cin_align_ind};
    cout_frates = d.frate{cout_align_ind};
    go_r        = d.trials.rat_dir==1;
    go_l        = d.trials.rat_dir==-1;
    psth_r           = cout_frates(go_r,:);
    psth_l           = cout_frates(go_l,:);
    
    nt = size(psth_r,2);
    
    parfor tt = 1:nt
        [cout_auc(cc,tt), cout_p(cc,tt), cout_ci(cc,tt,:)] = bootroc(psth_r(:,tt),psth_l(:,tt),nboot);
    end
end
%%
figure; 
imagesc(cout_auc,'x',cout_t)
caxis([0 1]+[1 -1].*.3)
line([0 0],ylim,'color','k')
colormap(colormapRedBlue)
colorbar
xlim([-1.75 .75])
%%
ax = example_cell_psth('cells',cellids(cc))
axpos = get(ax(2),'position')

ax(2).Position = axpos + [0 0 0 -.2]

ax2 = axes('position', [axpos(1) .8 axpos(3) .15])
plot(ax2, cout_t, cout_auc)
hold(ax2,'on')
plot(cout_t(cout_p < .05), .7 ,'k.')
line(xlim,[0 0]+.5)

xlim(ax2,get(ax(2),'xlim'))

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
