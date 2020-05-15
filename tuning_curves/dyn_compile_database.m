% [x, frbins, Pjoints, fr_given_as, fr_var_given_as, a_given_frs] = compile_database(cellid, t0s, {'force_dv', 0}, ...
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
%  force_dv     Default 0. If passed as 1, forces recomputation and resaving
%            in the database, even if the run existed before.
%
%  datadir   Default '~/Papers/TimHanks/PBupsPhys/Code/Carlosbin/compile_database_data'
%
%  alignment string matching one of align_strs in align_LUT. This will set
%            alignment of time zero of analysis (e.g., stimstart, cpokeend,
%            cpokeout)
%


function [x, frbins, Pjoints, fr_given_as, fr_var_given_as, a_given_frs] = dyn_compile_database(cellid, t0s,p, varargin)

frbins=0;
pairs = { ...
	'lag'        0.2         ; ...
	'frbins'     0:0.1:50    ; ...
	'dt'         0.01        ; ...
    'alignment'  'stimstart-cout-mask' ; ...      % string matching one of align_strs in align_LUT
	'direction'  'backward' ; ...       % 'forward' or 'backward'
	'trialnums'  []          ; ...
	'use_nans'   0           ; ...
	'datadir'    '/home/alex/Dropbox/spikes/cell_packager_data/compile_database_data'
	'force_dv'      0  ; ...
    'krn_width'  []          ; ...      % forces neural data to have this kernel width; if empty, uses whatever is in data    'krn_type'   []          ; ...      % forces neural data to have this kernel type; if empty, uses whatever is in data
    'krn_type'   'halfgauss' ; ...      % forces neural data to have this kernel type
    'fr_dt'      []          ; ...      % forces neural data to have this bin size; if empty, uses whatever is in data
    'norm_type'     'div'    ; ...      % type of firing rate normalization; 'div' is divisive, 'z' is z-score
    'fit_version'   'byrat'   ; ...      % fit params to use; 'byrat', 'combined', 'clickdiff'
    'n_iter'    1             ; ...      % number of iterations for refining estimate of DV
    'param_scale_num'        1     ; ... % parameter number to scale
    'param_scale_factor'     1     ; ... % multiplicative factor of that parameter
}; parseargs(varargin, pairs);


if ~exist(datadir, 'dir')
	mkdir(datadir);
end;
if ~exist([datadir filesep 'database.mat'], 'file')	
    
	index_columns = {'cellid' 't0s' 'lag' 'frbins' ...
		'dt' 'trialnums' 'use_nans' 'align_ind' ...
        'krn_width' 'fr_dt' 'direction' 'krn_type' ...
        'norm_type' 'fit_version' 'n_iterations' ...
        'param_scale_num' 'param_scale_factor'};
	
	database = cell(0, numel(index_columns));

	save([datadir filesep 'database'], 'index_columns', 'database');
end;

load([datadir filesep 'database']);

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

switch fit_version
    case 'byrat'
        fit_ind = 1;
    otherwise
        error('fit_version not recognized')
end

% index for database
% my_index = {cellid t0s lag frbins dt trialnums use_nans align_ind krn_width fr_dt dir_ind kt_ind nt_ind};
% my_index = {cellid t0s lag frbins dt trialnums use_nans align_ind krn_width fr_dt dir_ind kt_ind nt_ind fit_ind};
my_index = {cellid t0s lag frbins dt trialnums use_nans align_ind krn_width fr_dt dir_ind kt_ind nt_ind fit_ind n_iter param_scale_num param_scale_factor};

% Check if database is missing any columns from my_index
n_missing_cols = numel(my_index) - size(database,2);
if n_missing_cols > 0
    % if not, pad with default index values;
    % this assumes shorter databases use default value of 1 for the info specified by the newer columns; 
    % this is correct for the added columns so far (kt_ind, nt_ind, ft_ind, and n_iter);
    % would need to change code if different default were 
    database = cat(2,database, num2cell(ones(size(database,1),n_missing_cols)));
    
    % FIX: need to add param scaling factors here
	index_columns = {'cellid' 't0s' 'lag' 'frbins' ...
		'dt' 'trialnums' 'use_nans' 'align_ind' ...
        'krn_width' 'fr_dt' 'direction' 'krn_type' ...
        'norm_type' 'fit_version' 'n_iterations' ...
        'param_scale_num' 'param_scale_factor'};
end

u = find(cell_row_eq(database, my_index));
if force_dv || isempty(u) || ~exist([datadir filesep 'compiled_' num2str(u) '.mat'], 'file'),
    disp('Recompiling Joint distribution')
    % FIX: need to pass param scaling factor here
	[x, frbins, Pjoints, fr_given_as, fr_var_given_as, a_given_frs] = dyn_compile_dv2(cellid, t0s, p,'lag', lag, ...
		'frbins', frbins, 'dt', dt, 'trialnums', trialnums, 'use_nans', use_nans, 'alignment', alignment,...
        'krn_width', krn_width, 'fr_dt', fr_dt, 'krn_type', krn_type, 'norm_type', norm_type, ...
        'fit_version',fit_version,'n_iter',n_iter, 'param_scale_num', param_scale_num, ...
        'param_scale_factor', param_scale_factor);
	
	if isempty(u),
		u = size(database,1)+1;
		database = [database ; my_index]; %#ok<*NASGU>
	end;
    
    % make directory if it doesn't already exist
    if ~exist(datadir, 'dir'), mkdir(datadir); end;
    
	save([datadir filesep 'compiled_' num2str(u)], 'x', 'frbins', 'Pjoints', ...
		'fr_given_as', 'fr_var_given_as', 'a_given_frs', 'my_index');
	save([datadir filesep 'database'], 'index_columns', 'database');
else
	load([datadir filesep 'compiled_' num2str(u)]);
	return;
end;
