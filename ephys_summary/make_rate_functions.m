function [y, x, opts] = make_rate_functions(cellid, varargin)
% function [y, x, opts] = make_rate_functions(cellid, varargin)
%
% INPUT:
%
%   cellid: id of cell
%
% OUTPUT:
% 
% y                  trial X time matrix of smoothed rate functions
% x                  time vector for psth
%
% OVERRIDEABLE DEFAULTS:
%
% ex_cbreak_end = 1;      % excludes trials with NIC ending with legal_cbreak
% pre = 1;                % time before aligment event to include
% post = 1;               % time after aligment event to include
% krn_width = 0.1;        % time after aligment event to include
% bin_size = 0.025;       % time after aligment event to include
% % must be a string for a time event from vec_data
% % these include: 'cpoke_start', 'cpoke_end', 'cpoke_out', 'spoke_in', 'stim_start'
% ref_event = 'cpoke_end';
% post_mask_event = '';   % same strings as above
% pre_mask_event = '';    % same strings as above
% vec_data = [];          % pass in vec_data (must be passed array_data!)
% array_data = [];        % pass in array_data (must be passed with vec_data!)

% TO ADD:
%   1) other smoothing kernels




%% Default Parameters
p = inputParser;
addParameter(p, 'pre', 1)
addParameter(p, 'post', 1)
addParameter(p, 'krn_width', .1)
addParameter(p, 'krn_type', 'halfgauss')
addParameter(p, 'bin_size', .025)
addParameter(p, 'ref_event', 'cpoke_end')
addParameter(p, 'post_mask_event', '')
addParameter(p, 'pre_mask_event', '')
addParameter(p, 'pre_mask_offset', 0)
addParameter(p, 'post_mask_offset', 0)
addParameter(p, 'vec_data', [])
addParameter(p, 'array_data', [])
addParameter(p, 'mask_other_trials', [])
parse(p,varargin{:})
p = p.Results;
pre             = p.pre;                % time before aligment event to include
post            = p.post;               % time after aligment event to include
krn_width       = p.krn_width;        % kernel std in secs
krn_type        = p.krn_type; % either 'halfgauss' causal filter or 'fullgauss' acausal filter
bin_size        = p.bin_size;       % sampling distance after smoothing in secs
% must be a string for a time event from vec_data
% these include: 'cpoke_start', 'cpoke_end', 'cpoke_out', 'spoke_in', 'stim_start'
ref_event           = p.ref_event; 
post_mask_event     = p.post_mask_event;   % same strings as above
pre_mask_event      = p.pre_mask_event;    % same strings as above
post_mask_offset    = p.post_mask_offset;   % offset added to post-mask event
pre_mask_offset     = p.pre_mask_offset;    % offset added to pre-mask event
vec_data            = p.vec_data;          % pass in vec_data (must be passed array_data!)
array_data          = p.array_data;        % pass in array_data (must be passed with vec_data!)
mask_other_trials   = p.mask_other_trials;

% get the packaged phys data if it isn't passed already
if isempty(vec_data) || isempty(array_data)
    [array_data, vec_data] = package_pbups_phys(cellid);
end

% get reference times for alignment 
ref = eval(['vec_data.' ref_event]);
trial_start = -ref;
trial_dur   = vec_data.state_0_entries - vec_data.state_0_exits;
trial_end   = trial_start + trial_dur;

switch krn_type
    case 'halfgauss'
        % make half gaussian causal kernel
        dx=ceil(5*krn_width/bin_size);
        krn=normpdf(-dx:dx,0,krn_width/bin_size);
        krn(1:dx)=0;
        krn=(krn)/sum(krn)/bin_size;
        normalize_krn = 1;
    case 'fullgauss'
        % make half gaussian causal kernel
        dx=ceil(5*krn_width/bin_size);
        krn=normpdf(-dx:dx,0,krn_width/bin_size);
        krn=(krn)/sum(krn)/bin_size;
        normalize_krn = 1;
    case 'raster'
        krn = 1;
        normalize_krn = 0;
    otherwise
        error('Do not recognize krn_type')
end
        

% filter and align spikes
n_bins = ceil((pre+post)/bin_size);
y = nan*ones(length(array_data),n_bins);
for j=1:length(array_data)
    [y(j,:), x] = spike_filter(ref(j), array_data(j).spikes, krn,...
        'pre', pre, 'post', post, 'kernel_bin_size', bin_size,...
        'normalize_krn', normalize_krn);
end


% mask spikes; mask should be relative to ref because we are operating on
% output of spike_filter, which is ref-aligned
if ~isempty(post_mask_event) && ~isempty(pre_mask_event)
    post_mask   = eval(['vec_data.' post_mask_event]) - ref + post_mask_offset;
    pre_mask    = eval(['vec_data.' pre_mask_event]) - ref + pre_mask_offset;
    [y x]       = maskraster(x,y,pre_mask,post_mask);
elseif ~isempty(post_mask_event)
    post_mask   = eval(['vec_data.' post_mask_event]) - ref + post_mask_offset;
    [y x]       = maskraster(x,y,-Inf,post_mask);
elseif ~isempty(pre_mask_event)
    pre_mask    = eval(['vec_data.' pre_mask_event]) - ref + pre_mask_offset;
    [y x]       = maskraster(x,y,pre_mask,Inf);
end


if mask_other_trials
    [y x] = maskraster(x,y,-ref,trial_end);
    
end

