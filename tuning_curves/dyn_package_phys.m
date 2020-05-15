function [array_data vec_data sessid ratname] = dyn_package_phys(cellid, varargin)
% function [array_data vec_data sessid ratname] = package_pbups_phys(cellid, varargin)
%
% INPUT:
%
%   cellid: id of cell
%
%
% OUTPUT:
% 
% array_data will contain a struct array with the following fields, where each
% element of the array represents an individual, non-violation trial:
%   spikes:         times of spikes relative to state_0
%	left_bups:      times of bups played on left relative to state_0
%   right_bups:     times of bups played on right relative to state_0
%   head_dtheta:    [n x 2] matrix with columns timestamps and angular head velocity
%   parsed_events:  parsed events for the trial
%
%
% vec_data will contain a structure with fields that are n trial vectors with 
% the following elements:
%   cpoke_start:    time relative to state_0 when center poke is initiated
%   cpoke_end:      time relative to state_0 when cpoke1 state ends
%   cpoke_out:      time relative to state_0 when rat leaves center poke
%   spoke_in:       time relative to state_0 when rat enters side poke
%   stim_start:     time relative to state_0 of first click on either side
%   gamma:          generative gamma
%   pokedR:         true for right poke choice
%   state_0_exits:  state_0 exit times based on state machine clock
%   state_0_entries:  state_0 entry times based on state machine clock
%   bup_diff:       right bups minus left bups
%   good:           index to good (included) trials (e.g., no cpoke vio)
%
% sessid:           sessid for this cellid
% ratname:          name of rat for this cellid
%
% %% Default Parameters
% ex_cbreak_end = 1;      % excludes trials with NIC ending with legal_cbreak
% save_to_file = 0;       % saves packaged data to a file
% force = false;          % forces re-calculation even if file exists
% save_path = '~/ratter/svn_papers/TimHanks/PBupsPhys/PackagedData/'; % default save folder

% To Add: 
%   1) reward time


%% Default Parameters
p = inputParser();
addParamValue(p,'save_to_file',false);
addParamValue(p,'save_path','');
parse(p, varargin{:});


ex_cbreak_end = 1;      % excludes trials with NIC ending with legal_cbreak
force = false;          % forces re-calculation even if file exists

save_path = p.Results.save_path;
save_to_file = ~isempty(save_path);       % saves packaged data to a file

% override based on varargin
opts = overridedefaults(who, varargin);

% checks if data is already saved to file and loads it in that case
filename = [save_path 'dyn_phys_data_' num2str(cellid) '.mat'];
if ~force && exist(filename,'file')==2
    load(filename);
    return;
end

% get ephys data
[sessid, spikes] = bdata('select sessid, ts from spktimes where cellid="{S}"',cellid);
if isempty(sessid)
    fprintf(2,'Cell %d was not found',cellid);
    return;
end
spikes = spikes{1};

% get relevant data structures 
S = get_sessdata(sessid);
ratname = S.ratname{1};
pd = S.pd{1};
peh = S.peh{1};
protocol = S.protocol{1};
[pd,peh] = fix_sizes_in_pd(pd,peh);

% get state_0 times exit times
state_0_exits = extract_event(peh,'state_0(1,2)');
state_0_entries = extract_event(peh,'state_0(2,1)');

% get behavioral variables for each trial
% times relative to state_0
cpoke_start = extract_event(peh,'cpoke1(end,1)') - state_0_exits;
cpoke_end = extract_event(peh,'cpoke1(end,2)') - state_0_exits;
sstr = 'u = find(pk.C(:,1)<ps.wait_for_spoke(end,1),1,''last''); this_trial=pk.C(u,2)';
cpoke_out = extract_alignment(sessid, sstr, peh) - state_0_exits;
pokedR = (pd.sides=='r' & pd.hits==1) | (pd.sides=='l' & pd.hits==0);
spoke_in = nan * ones(size(cpoke_start));
gamma = nan * ones(size(cpoke_start));
bad_spoke = [];
for j = 1:length(peh),
    if ~isnan(pd.hits(j)) && isfield(peh(j).states, 'cpoke1') && ~isempty(peh(j).states.cpoke1)
        if pokedR(j)
            s_pokes = peh(j).pokes.R(:,1) - state_0_exits(j);
        else
            s_pokes = peh(j).pokes.L(:,1) - state_0_exits(j);
		end
		potential_side_poke = s_pokes(find(s_pokes > cpoke_end(j),1,'first'));
%  		potential_side_poke = s_pokes(find(s_pokes > cpoke_out(j),1,'first'));
		if isempty(potential_side_poke)
			% Side poke happened before the end of cpoke
			potential_side_poke = s_pokes(end);
            bad_spoke = [bad_spoke j];
		end;
		spoke_in(j) = potential_side_poke;
    end
    gamma(j) = pd.bupsdata{j}.gamma;
end;

if strcmpi(protocol, 'pbupsfc') | strcmpi(protocol, 'pbupsdyn'),
    for j = 1:length(peh),
        if ~isempty(peh(j).waves.sound_on)
            stim_start_delay(j) = peh(j).waves.sound_on(end,1) - state_0_exits(j) - cpoke_start(j);
        else
            stim_start_delay(j) = nan;
        end
    end
    stim_start_delay = stim_start_delay';
    good = find(pd.violations == 0);
    T = pd.samples;
else
    stim_start_delay = dyn_fetch_sph('StimulusSection_stim_start_delay', sessid);
    if strcmpi(protocol, 'samedifferent'),
        nvio = parse_cpoke_violations(peh);
        good = find(nvio == 0);
        
        T = nan * zeros(length(peh),1);
        for j = 1:length(peh),
            if ismember(j, good) && isfield(peh(j).states, 'cpoke1') && ~isempty(peh(j).states.cpoke1),
                T(j) = diff(peh(j).states.cpoke1(1,:)) - stim_start_delay(j);
            end;
        end;
    elseif strcmpi(protocol, 'pbups') ,
        good = find(pd.violations == 0);
        T = pd.samples;
    end;
end

% remove bad_spoke trials
try
good = setdiff(good, bad_spoke);

% exclude trials where nose is not in center when cpoke_timer pings
% this will exclude trials where NIC ends with legal_cbreak
if ex_cbreak_end
    if strcmpi(protocol, 'pbupsfc') | strcmpi(protocol, 'pbupsdyn'),
        for j=1:length(state_0_exits)
            if ~isempty(peh(j).waves.cpoke_timer)
                nic(j) = peh(j).waves.cpoke_timer(end,1)-cpoke_start(j)-state_0_exits(j);
            else
                nic(j) = nan;
            end
        end
        nic = nic';
    else
        nic = dyn_fetch_sph('StimulusSection_nose_in_center', sessid);
    end
    nic = nic(1:length(state_0_exits));
    timer_ping = cpoke_start + nic + state_0_exits;
    nic_ends_cin = nan * ones(size(cpoke_start));
    for j = 1:length(peh),
        c_pokes = peh(j).pokes.C;
        nic_ends_cin(j) = sum(c_pokes(:,1) < timer_ping(j) & c_pokes(:,2) > timer_ping(j))==1;
    end
    good = intersect(good, find(nic_ends_cin));
end
catch
    warning('couldn''t get rid of bad_spokes... maybe the protocol table doesn''t update properly');
end
% get head tracking data
try
    [hp_ts hp_data] = get_tracking(sessid);
    hp_ts = hp_ts(1:end-1); 
    dtheta = headvelocity(hp_ts, hp_data.theta)';
catch
    fprintf(2, 'Problem with head tracking info for session %d\n', sessid);
    dtheta = [];
end;

% loop through good trials and get the array data
g = 0; % need this to remove trials with no clicks
good_clicks = [];
for gi = 1:length(good),
    gtrial = good(gi);
    
    % clicks data
    left_bups = pd.bupsdata{gtrial}.left';
    left_bups = left_bups(left_bups < T(gtrial));
    
    right_bups = pd.bupsdata{gtrial}.right';
    right_bups = right_bups(right_bups < T(gtrial));
    
    % get time of first click
    first_left = nan;
    first_right = nan;
    if ~isempty(left_bups),  first_left  = left_bups(1); end;
    if ~isempty(right_bups), first_right = right_bups(1); end;
    
    % remove trials where no click is played
    if isnan(first_left) && isnan(first_right)
        continue
    else
        g = g+1;
        good_clicks = [good_clicks gi];
        vec_data.stim_start(g) = min([first_left first_right])+stim_start_delay(gtrial)+cpoke_start(gtrial);
        vec_data.bup_diff(g) = length(right_bups) - length(left_bups);
    end
    
    array_data(g).left_bups = left_bups + stim_start_delay(gtrial) + cpoke_start(gtrial);
    array_data(g).right_bups = right_bups + stim_start_delay(gtrial) + cpoke_start(gtrial);
    
    % neural and head position data
    sp_inds = spikes >= state_0_exits(gtrial) & spikes < state_0_entries(gtrial);
%     hp_inds = hp_ts >= state_0_exits(gtrial) & hp_ts < state_0_entries(gtrial);
    
    array_data(g).spikes = spikes(sp_inds) - state_0_exits(gtrial);
%     array_data(g).head_dtheta = [colvec(hp_ts(hp_inds))-state_0_exits(gtrial) ...
% 		colvec(dtheta(hp_inds))];
    
    % parsed events history
    array_data(g).parsed_events = peh(gtrial);
end

% remove bad "no click" trials from good list
good = good(good_clicks);

% get the good vector data
vec_data.stim_start = vec_data.stim_start';
vec_data.bup_diff = vec_data.bup_diff';
vec_data.cpoke_start = cpoke_start(good);
vec_data.cpoke_end = cpoke_end(good);
vec_data.cpoke_out = cpoke_out(good);
vec_data.spoke_in = spoke_in(good);
vec_data.gamma = gamma(good);
vec_data.pokedR = pokedR(good);
vec_data.state_0_exits = state_0_exits(good);
vec_data.state_0_entries = state_0_entries(good);
vec_data.good = good;
vec_data.T = T(good);

if save_to_file==1
    if ~exist(save_path, 'dir')
        mkdir(save_path);
    end;
    save(filename,'array_data','vec_data','sessid','ratname');
end
