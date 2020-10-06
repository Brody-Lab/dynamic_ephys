close all
clear all


dp = set_dyn_path;


warning('off','MATLAB:interp1:NaNinY');
warning('off','MATLAB:divideByZero');

%% get good cells 
cell_list = dyn_cells_db('force',0);   % That has to be run once to create cell_list
region      = 'fof';
select_str  = ['strcmp(region,''' region ''') &  normmean>1'];
select_str  = ['strcmp(region,''' region ''') '];

cellids     = cell2mat(extracting(cell_list, 'cellid', select_str));
rats     = cell2mat(extracting(cell_list, 'ratname', select_str));
alignment   = 'stimstart-cout-mask'; %'cpokeout'
lag         = 0;            %0.2;
t0s         = 0:0.025:1.9-lag;  % triggers rebinning and recompiling
n_dv_bins   = 100;          % triggers rebinning, but not recompiling 
krn_width   = 0.1;
force_frdv  = 1;            % keep as one if rebinning
force_bin   = 0;
force_dv    = 1;
norm_type   = 'none';
p.mask_after_stim_firing    = 1;
p.mask_stim_firing          = 0;
p.average_a_vals            = 1;
p.average_fr                = 1;
p.min_fr_mod                = 1;
p.fit_time_sigmoid          = 0;
p.use_nresidual             = 1;
p.plot_res                  = 3;
%cellid = [18181 18185 17799 17784 16875 18172 17870 16847 17855 17803 18466]
%cellid = [-1 -2 -3 -4 -5];
%cellids = [-1 -11 -2 -20 -21 -22 -23 -24 -3 -30 -31 -32];
results = dyn_fr_dv_map(cellids, t0s, n_dv_bins, p,...
    'lag', lag, 'krn_width', krn_width, 'alignment', alignment, ...
    'var_weight', false, 'force_frdv',force_frdv,'force_bin',force_bin,...
    'force_dv', force_dv,'norm_type',norm_type);

%%
if 0
good_cells = results.cell_dex;
good_rats = rats(good_cells,:);
counts = zeros(4,2);
counts(1,1) = sum(all(ismember(rats,'H037'),2));
counts(2,1) = sum(all(ismember(rats,'H066'),2));
counts(3,1) = sum(all(ismember(rats,'H084'),2));
counts(4,1) = sum(all(ismember(rats,'H129'),2));
counts(1,2) = sum(all(ismember(good_rats,'H037'),2));
counts(2,2) = sum(all(ismember(good_rats,'H066'),2));
counts(3,2) = sum(all(ismember(good_rats,'H084'),2));
counts(4,2) = sum(all(ismember(good_rats,'H129'),2));
h037 = all(ismember(rats,'H037'),2);
h066 = all(ismember(rats,'H066'),2);
h084 = all(ismember(rats,'H084'),2);
h129 = all(ismember(rats,'H129'),2);

results = compute_variance_explained(results, 18185,1)
figure; hold on;
%celldex = find(results.cellid == 18181);
%celldex = 2;
%plot_cell(results, results.cellid(celldex));
%figure; plot(results.t0s, results.frm_time(:,celldex), 'k');
%figure; plot(results.t0s, results.var_explained(:,celldex),'r')

results = decompose_tuning_curves(results);
plot_decomposition_analysis(results)
end


