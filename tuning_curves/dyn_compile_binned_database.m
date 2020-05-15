% [x, frbins, Pjoints, fr_given_as,  fr_var_given_as, a_given_frs] = compile_binned_database(cellid, t0s, n_dv_bins, {'force', 0}, ...
%         {'datadir', '~/svn_papers/TimHanks/PBupsPhys/Code/Carlosbin/compile_database_data'}, ...
%         {'lag', 0.2}, {'frbins', 0:0.1:10}, {'dt', 0.01}, {'trialnums', []}))
%
% Maintains a database of runs of compile_dv.m using the indicated cellid,
% t0s, and optional params. If a run in the database exists, loads that;
% otherwise runs compile_dv.m to compute the answer (and save it in the
% database).
%
% See compile_dv.m for parameters, optional parameters, and return
% parameters.
%
% In addition, there are these two optional parameters:
%
%  force     Default 0. If passed as 1, forces recomputation and resaving
%            in the database, even if the run existed before.
%
%  datadir   Default '~/svn_papers/TimHanks/PBupsPhys/Code/Carlosbin/compile_binned_database_data'
%
%  alignment string matching one of align_strs in align_LUT. This will set
%            alignment of time zero of analysis (e.g., stimstart, cpokeend,
%            cpokeout)
%



function [x, frbins, Pjoints, fr_given_as, fr_var_given_as, a_given_frs] = dyn_compile_binned_database(cellid, t0s, n_dv_bins, p,varargin)

frbins=0;
pairs = { ...
	'lag'        0.2         ; ...
	'frbins'     0:0.1:50    ; ...
	'dt'         0.01        ; ...
    'alignment'      'stimstart-cout-mask' ; ...      % string matching one of align_strs in align_LUT
	'direction'      'backward' ; ...       % 'forward' or 'backward'
	'trialnums'  []          ; ...
	'use_nans'   0           ; ...
	'datadir'          '~/Dropbox/spikes/cell_packager_data'
	'force_bin'              0 ; ...
    'force_dv'               0; ...
    'krn_width'  []          ; ...      % forces neural data to have this kernel width; if empty, uses whatever is in data    'krn_type'   []          ; ...      % forces neural data to have this kernel type; if empty, uses whatever is in data
    'krn_type'   'halfgauss' ; ...      % forces neural data to have this kernel type
    'fr_dt'      []          ; ...      % forces neural data to have this bin size; if empty, uses whatever is in data
    'norm_type'     'div'    ; ...      % type of firing rate normalization; 'div' is divisive, 'z' is z-score
    'fit_version'   'byrat'   ; ...      % fit params to use; 'byrat', 'combined', 'clickdiff'
    'n_iter'    1            ; ...      % number of iterations for refining estimate of DV
    'param_scale_num'        1     ; ... % parameter number to scale
    'param_scale_factor'     1     ; ... % multiplicative factor of that parameter
}; parseargs(varargin, pairs);


% error checking
if n_dv_bins<3
    error('You must have at least three DV bins. 1 for each bound and 1 for non-bound values.');
end

if ~exist(datadir, 'dir')
	mkdir(datadir);
end;
if ~exist([datadir filesep 'database.mat'], 'file')	
	index_columns = {'cellid' 't0s' 'lag' 'frbins' ...
		'dt' 'trialnums' 'use_nans' 'align_ind' ...
        'krn_width' 'fr_dt' 'n_dv_bins' 'direction' ...
        'krn_type' 'norm_type' 'fit_version' 'n_iterations' ...
        'param_scale_num' 'param_scale_factor'};
	
	database = cell(0, numel(index_columns));

	save([datadir filesep 'database'], 'index_columns', 'database');
end;

load([datadir filesep 'database']);

Pjoints     = [];
fr_given_as = [];
fr_var_given_as = [];

% find alignment index
align_strs = align_LUT;
align_ind = strmatch(alignment,align_strs,'exact');


% make direction into index for databasing
switch direction
    case 'forward'
        dir_ind = 1;
    case 'backward'
        dir_ind = 2;
    otherwise
        error('Direction not recognized')
end

% make krn_type into index for databasing
switch krn_type
    case 'halfgauss'
        kt_ind = 1;
    case 'fullgauss'
        kt_ind = 2;
    otherwise
        error('krn_type not recognized')
end

% make norm_type into index for databasing
switch norm_type
    case 'div'
        nt_ind = 1;
    case 'z'
        nt_ind = 2;
    case 'none'
        nt_ind = 3;
    otherwise
        error('norm_type not recognized')
end

% make fit version into index for databasing
switch fit_version
    case 'byrat'
        fit_ind = 1;
    case 'combined'
        fit_ind = 2;
    case 'byrat2'
        fit_ind = 3;
    case 'deltaclicks'
        fit_ind = 4;
    case 'byrat3'
        fit_ind = 5;
    case 'byrat35'
        fit_ind = 35;
    case 'byrat35v2'
        fit_ind = 36;
    case 'combined35'
        fit_ind = 37;
    case 'byrat35_ncd'
        fit_ind = 38;
    case 'byrat35_phys'
        fit_ind = 39;
    otherwise
        error('fit_version not recognized')
end

% index for database
% my_index = {cellid t0s lag frbins dt trialnums use_nans align_ind krn_width fr_dt n_dv_bins dir_ind kt_ind nt_ind fit_ind};
my_index = {cellid t0s lag frbins dt trialnums use_nans align_ind krn_width fr_dt n_dv_bins dir_ind kt_ind nt_ind fit_ind n_iter param_scale_num param_scale_factor};

% Check if database is missing any columns from my_index
n_missing_cols = numel(my_index) - size(database,2);
if n_missing_cols > 0
    % if not, pad with default index values;
    % this assumes shorter databases use default value of 1 for the info specified by the newer columns; 
    % this is correct for the added columns so far (kt_ind and nt_ind);
    % would need to change code if different default were 
    database = cat(2,database, num2cell(ones(size(database,1),n_missing_cols)));
    
	index_columns = {'cellid' 't0s' 'lag' 'frbins' ...
		'dt' 'trialnums' 'use_nans' 'align_ind' ...
        'krn_width' 'fr_dt' 'n_dv_bins' 'direction' ...
        'krn_type' 'norm_type' 'fit_version' 'n_iterations' ...
        'param_scale_num' 'param_scale_factor'};
end

u = find(cell_row_eq(database, my_index));
if force_bin || isempty(u) || ~exist([datadir filesep 'compiled_binned_' num2str(u) '.mat'], 'file'),
    disp('Rebinning')
	[x, frbins, Pjoints_unbinned, ~, ~, a_given_frs] = dyn_compile_database(cellid, t0s, p,'lag', lag, ...
		'frbins', frbins, 'dt', dt, 'trialnums', trialnums, 'use_nans', use_nans, 'alignment', alignment,...
        'krn_width', krn_width, 'fr_dt', fr_dt, 'force_dv', force_dv, 'direction', direction, ...
        'krn_type', krn_type, 'norm_type', norm_type,'fit_version',fit_version,'n_iter',n_iter, ...
        'param_scale_num', param_scale_num, 'param_scale_factor', param_scale_factor);

    % calculate binned 'Pjoints', 'fr_given_as', 'fr_var_given_as' for the given n_dv_bins
    x_bound_margin = (x(2)-x(1)) * 0.5;
    dv_bin_edges = linspace(min(x)+x_bound_margin,max(x)-x_bound_margin,n_dv_bins-1);
    dv_axis = dv_bin_edges + 0.5*(dv_bin_edges(2)-dv_bin_edges(1));
    dv_axis = dv_axis(1:end-1);
    dv_axis = [min(x) dv_axis max(x)];  % add bound bins for axis
    
    % FIX: account for behavioral bias
    %     dv_bin_edges = dv_bin_edges + rat_bias;

    if numel(dv_axis) > size(Pjoints_unbinned,3)
        error('You are binning accumulation values more finely than the model generates')
    end
    for time_i=1:numel(t0s)
        if isempty(Pjoints),
            if use_nans,
                Pjoints     = NaN*ones(numel(t0s), numel(frbins), numel(dv_axis));
            else
                Pjoints     = zeros(numel(t0s), numel(frbins), numel(dv_axis));
            end;
            fr_given_as = zeros(numel(t0s), numel(dv_axis));
            fr_var_given_as = zeros(numel(t0s), numel(dv_axis));
        end;

        % make binned Pjoints
        % Deal with the first bin
        these_dv_bins       = x < dv_bin_edges(1);
        Pjoints(time_i,:,1) = nansum(Pjoints_unbinned(time_i,:,these_dv_bins),3);
        edge_match          = x == dv_bin_edges(1);
        Pjoints(time_i,:,1) = Pjoints(time_i,:,1) + 0.5*nansum(Pjoints_unbinned(time_i,:,edge_match),3); 

        % Deal with middle bins
        for j=2:n_dv_bins-1
            these_dv_bins       = x>dv_bin_edges(j-1) & x<dv_bin_edges(j);
            Pjoints(time_i,:,j) = nansum(Pjoints_unbinned(time_i,:,these_dv_bins),3); % Pjoints is t0 x FR x dv
            % deal with bin edges that exactly match DV centers; for
            % these case, split joint dist in half between adjacent bins
            edge_match          = x==dv_bin_edges(j-1) | x==dv_bin_edges(j);
            Pjoints(time_i,:,j) = Pjoints(time_i,:,j) + 0.5*nansum(Pjoints_unbinned(time_i,:,edge_match),3); 
        end
        % Deal with last bin       
        these_dv_bins           = x > dv_bin_edges(end);
        Pjoints(time_i,:,end)   = nansum(Pjoints_unbinned(time_i,:,these_dv_bins),3);
        edge_match              = x == dv_bin_edges(end);
        Pjoints(time_i,:,end)   = Pjoints(time_i,:,end) + 0.5*nansum(Pjoints_unbinned(time_i,:,edge_match),3); 

        if abs(sum(sum(Pjoints(time_i,:,:))) - sum(sum(Pjoints_unbinned(time_i,:,:)))) > 1e-6
            keyboard
            error('Probability mass is leaking during binning')
        end
 
        myPj = squeeze(Pjoints(time_i,:,:)); 
        Pjoint_given_a = myPj ./ (ones(size(myPj,1),1)*sum(myPj,1));
        fr_given_as(time_i,:) = Pjoint_given_a'*frbins';
        fr_var_given_as(time_i,:) = (Pjoint_given_a'*(frbins'.^2)) - (fr_given_as(time_i,:)'.^2);
    end;

	if isempty(u),
		u = size(database,1)+1;
		database = [database ; my_index]; %#ok<*NASGU>
	end;
    
    % make directory if it doesn't already exist
    if ~exist(datadir, 'dir'), mkdir(datadir); end;
    
    % save binned dv_axis as x
    x = dv_axis;
    
	save([datadir filesep 'compiled_binned_' num2str(u)], 'x', 'frbins', 'Pjoints', ...
		'fr_given_as', 'fr_var_given_as', 'a_given_frs', 'my_index');
	save([datadir filesep 'database'], 'index_columns', 'database');
else
	load([datadir filesep 'compiled_binned_' num2str(u)]);
	return;
end;

