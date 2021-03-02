function [res,cellids,computed, not_computed] = ...
    load_STA_population(which_switch, varargin)

p = inputParser;
addParameter(p,'force', 0)
addParameter(p,'slim_data',1)
addParameter(p,'bad_strength',0)
addParameter(p,'min_t',-Inf)
addParameter(p,'max_t',Inf)
parse(p,varargin{:})

p = p.Results;
this_force = p.force;
slim_data = p.slim_data;
bad_strength = p.bad_strength;
min_t = p.min_t;
max_t = p.max_t;

select_str = 'normmean > 0';
dp = set_dyn_path;
sta_fn = ['population_STA_' which_switch '_' num2str(bad_strength) '.mat'];
sta_file = fullfile(dp.spikes_dir, sta_fn);
if ~this_force && exist(sta_file,'file')
    load(sta_file);
else
    % get the list of cells
    cell_list = dyn_cells_db;
    % only look at cells with reasonable firing rates during the trial
    
    save_switch = [which_switch '_' num2str(bad_strength)];

    %%%%% logic here for rat specific model
    if which_switch(end-4) == '_' % Rat specific
        this_rat = which_switch(end-3:end);
        which_switch = which_switch(1:end-5);
        select_str = ['normmean > .5 & strcmp(ratname, ''' this_rat ''')'];
    end 

    cellids = cell2mat(extracting(cell_list, 'cellid',select_str));
%    sessids = cell2mat(extracting(cell_list, 'sessid',select_str));
    % iterate over cells and compute STAs
    computed        = zeros(size(cellids));
    not_computed    = ones(size(cellids));
    res = cell(1,length(cellids));
    nn = 1;
     for cc = nn:length(cellids)
         try
             res{nn} =  compute_switch_triggered_average(cellids(cc),...
                 'post',2,'which_switch',which_switch, ...
                 'n_shuffles', 1000,'save_file',1,...
                 'mask_other_switch',1, 'bad_strength', bad_strength);
             res{nn}.STR_right_shuff = [];
             res{nn}.STR_left_shuff  = [];
             res{nn}.STR_right_real  = nanmean(res{nn}.STR_right_real);
             res{nn}.STR_left_real   = nanmean(res{nn}.STR_left_real);
             res{nn}.dprime_shuff    = [];
            if slim_data
             res{nn}.params          = [];
            end
             disp(cc)
             nn = nn+1;
            computed(cc) = 1;
            not_computed(cc) = 0;
         catch ME
            disp(ME.message)
         end 
     end
     
    save(sta_file,'cellids','cc','res','not_computed','computed')
end


