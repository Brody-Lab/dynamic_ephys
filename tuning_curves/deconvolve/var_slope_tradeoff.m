close all;
clear all;
cd ~/Dropbox/spikes/bin/tuning_curves/

addpath ~/ratter/svn_papers/TimHanks/PBupsPhys/Code/
addpath ~/ratter/svn_papers/TimHanks/PBupsPhys/Code/Carlosbin
addpath ~/ratter/ExperPort/bin
addpath ~/ratter/Analysis/Pbups
addpath ~/ratter/ExperPort/MySQLUtility
addpath ~/ratter/ExperPort/Analysis
addpath ~/ratter/ExperPort/Analysis/SameDifferent/
addpath ~/ratter/ExperPort/HandleParam
addpath ~/ratter/Analysis/helpers
addpath ~/Dropbox/spikes/cell_packager_data/
addpath ~/Dropbox/spikes/bin
addpath ~/Dropbox/spikes/bin/tuning_curves

warning('off','MATLAB:interp1:NaNinY');
warning('off','MATLAB:divideByZero');

%% get good cells 
cell_list = dyn_cells_db('force',0);   % That has to be run once to create cell_list
region      = 'fof';
select_str  = ['strcmp(region,''' region ''') &  normmean>1'];
cellids     = cell2mat(extracting(cell_list, 'cellid', select_str));

alignment   = 'stimstart-cout-mask'; %'cpokeout'
t0s         = 0:0.025:1.9;  % triggers rebinning and recompiling
n_dv_bins   = 100;          % triggers rebinning, but not recompiling 
lag         = 0;            %0.2;
krn_width   = 0.1;
force_frdv  = 1;            % keep as one if rebinning
force_bin   = 1;
force_dv    = 1; 
p.mask_after_stim_firing    = 1;
p.mask_stim_firing          = 0;
p.average_a_vals            = 1;
p.average_fr                = 1;

cellid = [-4 -6 -7 -8 -9 -10];
results     = dyn_fr_dv_map(cellid, t0s, n_dv_bins, p,'lag', lag, 'krn_width', krn_width, 'alignment', alignment, 'var_weight', false, 'force_frdv',force_frdv,'force_bin',force_bin,'force_dv', force_dv);




cd ~/Dropbox/spikes/bin/tuning_curves/deconvolve/
vars = zeros(size(cellid));
for i=1:length(cellid)
    load(['var_data_' num2str(abs(cellid(i))) '.mat'])
    vars(i) = average_var;
end

gen_slopes = zeros(size(cellid));
for i=1:length(cellid)
    x = results.dv_axis;
    big_x = -40:max(diff(x)):40;
    step = double(big_x > 0);
    filter = normpdf(big_x,0,vars(i));
    sig = conv(step,filter,'same');
    sig = sig - min(sig);
    sig = sig./max(sig);
    [betas,resid,jacob,sigma,mse] = nlinfit(big_x,sig,@dyn_sig,[0, 1/3]);
    gen_slopes(i) = betas(2)/4;
    [yhat,H] = wienerFilter(step,sig);
end


figure; plot(vars,results.slope_cell,'bo-','markerfacecolor','b')
hold on;

%plot(vars,gen_slopes, 'mo-','markerfacecolor','r')
plot(vars(1), results.slope_cell(1), 'ro-','markersize',20)
legend({'Tuning Curve','Convolved Step function','data variance'})
ylabel('Recovered Slope of sigmoid')
xlabel('Average model variance')
set(gca, 'fontsize',12)
xlim([0 max(vars)+1])
ylim([0 1])

