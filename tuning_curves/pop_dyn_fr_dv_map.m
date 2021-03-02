function res = pop_dyn_fr_dv_map(cellids, sessid, varargin)
% function res = dyn_fr_dv_map(cellid, varargin)
p = inputParser;
addParameter(p,'lag',        0);
addParameter(p,'frbins',    []);
addParameter(p,'use_nans',   0);
addParameter(p,'mask_after_stim',   true);
addParameter(p,'mask_during_stim',  false);
addParameter(p,'average_a_vals',    true)
addParameter(p,'average_fr',        true);
addParameter(p,'end_mask_s',        0);
addParameter(p,'dt',        0.01);
addParameter(p,'max_t',     1.9);
addParameter(p,'t0s',    []);
addParameter(p,'alignment', 'stimstart');   % string matching one of align_strs in dyn_align_LUT
addParameter(p,'trialnums',  []);
addParameter(p,'krn_width',  []);      % forces neural data to have this kernel width; if empty, uses whatever is in data    'krn_type'   []          );      % forces neural data to have this kernel type; if empty, uses whatever is in data
addParameter(p,'krn_type',   'halfgauss');      % forces neural data to have this kernel type
addParameter(p,'fr_dt',      0.0250);      % forces neural data to have this bin size; if empty, uses whatever is in data
addParameter(p,'norm_type',     'none');      % type of firing rate normalization; 'div' is divisive, 'z' is z-score, 'none' no normalization
addParameter(p,'frates',   []);
addParameter(p,'shuffle_trials',   0);
addParameter(p,'model',   []);
addParameter(p,'n_dv_bins',   []);
addParameter(p,'which_switch',   []);
addParameter(p,'demean_frates',   0);
addParameter(p,'zscore_frates',   0);
parse(p,varargin{:});
p = p.Results;


if isempty(sessid)
    sessid = zeros(size(cellids));
    for cc = 1:length(cellids)
        sessid(cc) = bdata('select sessid from cells where cellid={S}',cellids(cc));
    end
end

unique_sessid = unique(sessid);
fprintf('session (of %i) ...',length(unique_sessid))
for ss = 1:length(unique_sessid)
    fprintf('%i...',ss)
    this_sessid = unique_sessid(ss);
    ind     = find(sessid == unique_sessid(ss));
    model   = get_data_model_p(this_sessid);
    for ff = ind;
        res(ff) = dyn_fr_dv_map(cellids(ff), ...
            varargin{:}, 'model', model);
    end
end
    
