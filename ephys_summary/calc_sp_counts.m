function [counts, durations, opts] = calc_sp_counts(cellid, varargin)
% function [y, x, opts] = calc_sp_counts(cellid, varargin)
%
% INPUT:
%
%   cellid: id of cell
%
% OUTPUT:
% 
% counts            trial vector of spike counts within mask window
% durations         duration of mask window             	
%
% OVERRIDEABLE DEFAULTS:
%
% ex_cbreak_end = 1;      % excludes trials with NIC ending with legal_cbreak
% pre = 1;                % time before aligment event to include
% post = 1;               % time after aligment event to include
% % must be a string for a time event from vec_data
% % these include: 'cpoke_start', 'cpoke_end', 'cpoke_out', 'spoke_in', 'stim_start'
% ref_event = 'cpoke_end';
% post_mask_event = '';   % same strings as above
% pre_mask_event = '';    % same strings as above
% vec_data = [];          % pass in vec_data (must be passed array_data!)
% array_data = [];        % pass in array_data (must be passed with vec_data!)


%% Default Parameters
p = inputParser;
addParameter(p, 'ex_cbreak_end', 1)
addParameter(p, 'pre', 1)
addParameter(p, 'post', 1)
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
ex_cbreak_end   = p.ex_cbreak_end;      % excludes trials with NIC ending with legal_cbreak
pre             = p.pre;                % time before aligment event to include
post            = p.post;               % time after aligment event to include
% must be a string for a time event from vec_data
% these include: 'cpoke_start', 'cpoke_end', 'cpoke_out', 'spoke_in', 'stim_start'
ref_event           = p.ref_event; 
post_mask_event     = p.post_mask_event;   % same strings as above
pre_mask_event      = p.pre_mask_event;    % same strings as above
post_mask_offset    = p.post_mask_offset;   % offset added to post-mask event
pre_mask_offset     = p.pre_mask_offset;    % offset added to pre-mask event
vec_data            = p.vec_data;          % pass in vec_data (must be passed array_data!)
array_data          = p.array_data;        % pass in array_data (must be passed with vec_data!)

% get the packaged phys data if it isn't passed already
if isempty(vec_data) || isempty(array_data)
    [array_data, vec_data] = package_dyn_phys(cellid, ...
        'ex_cbreak_end', ex_cbreak_end);
end

% get reference times to calcuate pre,post times
ref = eval(['vec_data.' ref_event]);
ref_pre = ref - pre;
ref_post = ref + post;

% mask times
if isempty(pre_mask_event)
    % no pre-mask
    pre_mask = 0;   
else
    pre_mask = eval(['vec_data.' pre_mask_event]) + pre_mask_offset;
end
if isempty(post_mask_event)
    % no post-mask; limit to end of trial
    post_mask = vec_data.state_0_entries - vec_data.state_0_exits;
else
    post_mask = eval(['vec_data.' post_mask_event]) + post_mask_offset;
end

% start and end times 
start_time = max(ref_pre, pre_mask);
end_time = min(ref_post, post_mask);
durations = end_time - start_time;

counts = nan(size(durations));
for i=1:length(array_data)
    try
        counts(i) = qcount(array_data(i).spikes, start_time(i), end_time(i));
    catch
        counts(i) = 0;
    end
end
