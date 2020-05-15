function cell_list = dyn_cells_db(varargin)
% Function that makes databases of info on pbups cells by running
% cell_packager. The "database" is saved as structure called "cell_list".
% There is a row for each cell plus a top row that gives column headings.
% Each column saves one scalar value for each cell. These are used for cell
% selection for further analysis.


% Default Params 
force = 0; 
force_each = 0;
repack_each = 0;
save_dir = '~/Dropbox/spikes/cell_packager_data'; 
save_db = true;  
% override based on varargin
opts = overridedefaults(who, varargin);
save_name = 'dyn_multi_db.mat';

save_path = fullfile(save_dir, save_name);


if ~force && exist(save_path, 'file')
    load(save_path);
else
    
        
    [~, multi] = dyn_get_cells()
    
    cell_list = {'cellnum' 'region' 'ratname' 'cellid' 'sessid' 'sessiondate' 'lr' 'hitfrac' 'mingamma' 'brainsideright' ...
        'prefsideright' 'prefp' 'normmean' 'prefsideright_mv' 'prefp_mv' 'min_p' 'prefsideright_cout100' 'prefp_cout100' ...
        'prefsideright_stimcout' 'prefp_stimcout'...
        'prefoutcomehit' 'hitprefp'...
        'prefoutcomehit_mv' 'hitprefp_mv'...
        'prefoutcomehit_cout100' 'hitprefp_cout100'...
        'prefoutcomehit_stimcout' 'hitprefp_stimcout'};
    
    % load all single units
    %singles = get_cells('pbups_all');
    
    
    for i=1:numel(multi)
        fprintf('Packaging cell %d of %d \n', i, numel(multi));
%         data = cell_packager(singles(i),'krn_type','fullgauss','krn_width',0.05);

        data = dyn_cell_packager(multi(i),'krn_type','halfgauss','krn_width',0.1,...
            'datadir' , save_dir, 'force', force_each, 'repack', repack_each);
        stim_index = strmatch('stimend',data.align_strs,'exact');
        stim_cout_index = strmatch('stimstart-cout-mask',data.align_strs,'exact');
        mv_index = strmatch('postmove',data.align_strs,'exact');
        cout100_index = strmatch('cout100',data.align_strs,'exact');
        stats_index = strmatch('cpokeend-choice-all',data.stats_strs,'exact');
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
    
    if save_db,
        if ~exist(save_dir, 'dir'), mkdir(save_dir); end;
        save(save_path, 'cell_list');
    end
end


