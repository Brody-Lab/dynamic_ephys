function cell_list = dyn_cells_db(varargin)
% Function that makes databases of info on pbups cells by running
% cell_packager. The "database" is saved as structure called "cell_list".
% There is a row for each cell plus a top row that gives column headings.
% Each column saves one scalar value for each cell. These are used for cell
% selection for further analysis.

p = inputParser();
addParameter(p, 'force', false)
addParameter(p, 'repack_each', false)
addParameter(p, 'save_dir', [])
addParameter(p, 'save_name', [])
addParameter(p, 'dyn_path', [])
addParameter(p, 'do_save', true)
addParameter(p, 'krn_type', 'halfgauss');
addParameter(p, 'krn_width', .1);
parse(p, varargin{:});
par     = p.Results;

if isempty(par.dyn_path)
    dp      = set_dyn_path;
else
    dp      = par.dyn_path;
end

force       = par.force; 
repack_each = par.repack_each;
save_dir    = par.save_dir; 
save_name   = par.save_name; 
do_save     = par.do_save;  
krn_type    = par.krn_type;
krn_width   = par.krn_width;

if isempty(par.save_name)
    save_name = dp.celldat_filename;
end
if isempty(par.save_dir)
    save_dir = dp.celldat_dir;
end

save_path   = fullfile(save_dir, save_name);

if ~exist(save_dir,'dir')
    in = input(sprintf(['Save directory doesn''t exist\n'...
        '%s\n would you like to create it and continue (y/n)?'],save_dir), 's');
    if lower(in) == 'y'
        mkdir(save_dir);
    else
        in2 = input('continue without saving? (y/n)','s');
        if lower(in2) == 'y'
            warning('continuing without saving')
        else
            warning('quitting without retrieving cells');
            return
        end
    end  
end

if ~force && exist(save_path, 'file')
    warning('loading database from file. Will not update new cells')
    load(save_path);
else
    fprintf('Attempting to connect to bdata');
    bdata('connect');
    fprintf('Getting cell list\n');
    [singles, multi] = dyn_get_cells();
    
    cell_list = {'cellnum' 'region' 'ratname' 'cellid' 'sessid' 'sessiondate' 'lr' 'hitfrac' 'mingamma' 'brainsideright' ...
        'prefsideright' 'prefp' 'normmean' 'prefsideright_mv' 'prefp_mv' 'min_p' 'prefsideright_cout100' 'prefp_cout100' ...
        'prefsideright_stimcout' 'prefp_stimcout'...
        'prefoutcomehit' 'hitprefp'...
        'prefoutcomehit_mv' 'hitprefp_mv'...
        'prefoutcomehit_cout100' 'hitprefp_cout100'...
        'prefoutcomehit_stimcout' 'hitprefp_stimcout'};
    
    for i=1:numel(multi)
        fprintf('Packaging cell %d of %d \n', i, numel(multi));

        data = dyn_cell_packager(multi(i),'krn_type',krn_type,'krn_width',krn_width,...
            'datadir' , save_dir, 'repack', repack_each);
        stim_index      = strmatch('stimend',data.align_strs,'exact');
        stim_cout_index = strmatch('stimstart-cout-mask',data.align_strs,'exact');
        mv_index        = strmatch('postmove',data.align_strs,'exact');
        cout100_index   = strmatch('cout100',data.align_strs,'exact');
        stats_index     = strmatch('cpokeend-choice-all',data.stats_strs,'exact');
        cell_list = [cell_list ; ...
            {i data.region data.ratname data.cellid data.sessid data.sessiondate ...
            data.lapse_rate mean(data.trials.hit) data.min_gamma data.brainsideright ...
            data.prefsideright{stim_index} data.prefp{stim_index} data.norm_mean ...
            data.prefsideright{mv_index} data.prefp{mv_index} data.min_p{stats_index} ...
            data.prefsideright{cout100_index} data.prefp{cout100_index} ...
            data.prefsideright{stim_cout_index} data.prefp{stim_cout_index} ...
            data.prefoutcomehit{stim_index} data.hitprefp{stim_index}  ...
            data.prefoutcomehit{mv_index} data.hitprefp{mv_index}  ...
            data.prefoutcomehit{cout100_index} data.hitprefp{cout100_index} ...
            data.prefoutcomehit{stim_cout_index} data.hitprefp{stim_cout_index}}]; %#ok<*AGROW>
    end;
    
    if do_save
        if ~exist(save_dir, 'dir'), mkdir(save_dir); end;
        save(save_path, 'cell_list');
    end
end


