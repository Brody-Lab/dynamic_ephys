function [switch_to_0, switch_to_1, array_data, vec_data] = ...
    get_switches(cellid, varargin)

p = inputParser;
addParameter(p, 'which_switch', 'model')
addParameter(p, 'array_data', [])
addParameter(p,'vec_data',[]);
addParameter(p,'clear_bad_strengths',1);
addParameter(p,'bad_strength',0);
addParameter(p,'fit_line',1);
addParameter(p,'exclude_final',0);
addParameter(p,'final_only',0);
addParameter(p,'min_pre_dur',0);
addParameter(p,'min_post_dur',0);
addParameter(p,'min_switch_t',0);
addParameter(p,'max_switch_t',2);
addParameter(p,'which_trials',[]);
parse(p,varargin{:});
p = p.Results;

array_data  = p.array_data;
vec_data    = p.vec_data;

if isempty(array_data) || isempty(vec_data)
    [array_data, vec_data, this_sessid] = package_dyn_phys(cellid);
else
    this_sessid = bdata('select sessid from cells where cellid={S}',cellid);
    this_sessid = this_sessid{1};
end

% align events to stimulus onset
array_data = cleanup_array_data(array_data, vec_data);
% add generative state switches to array_data
array_data = compute_state_switches(array_data);
% decide either to use the model or generative switches
switch p.which_switch
    case 'model'
        % if array_data wasn't passed in with model switches, compute them
        % using the saved model mean and the relevant model parameters
        if ~isfield(array_data, 'model_switch_to_0') ||  ~isfield(array_data, 'model_switch_to_1')
            array_data = compute_model_switches(array_data, this_sessid, ...
                'remove_initial_choice', 1, 'eval_dt', 1e-3, ...
                'strength_window', .1, ...
                'clear_bad_strengths', p.clear_bad_strengths, ...
                'bad_strength', p.bad_strength, 'fit_line', p.fit_line);
        end
        switch_to_0 = {array_data.model_switch_to_0};
        switch_to_1 = {array_data.model_switch_to_1};
        %stim_start = vec_data.stim_start(includes);
    case 'generative'
        switch_to_0 = {array_data.switch_to_0};
        switch_to_1 = {array_data.switch_to_1};
        %stim_start = vec_data.stim_start(includes);
    case 'accumulation'
        array_data = compute_accumulation_state_switches(array_data);
        switch_to_0 = {array_data.accumulation_switch_to_0};
        switch_to_1 = {array_data.accumulation_switch_to_1};
    otherwise
        error('which_switch improperly specified')
end

which_trials = p.which_trials;
if isempty(which_trials)
    which_trials = true(size(vec_data.good));
end
assert(isequal(length(which_trials),length(switch_to_0)));
switch_to_0 = switch_to_0(which_trials);
switch_to_1 = switch_to_1(which_trials);

T = vec_data.stim_dur;
NT = length(switch_to_0);
for tt = 1:NT
   this_to_0 = sort(switch_to_0{tt}); 
   this_to_1 = sort(switch_to_1{tt});
   if ~isempty(this_to_0)
       isi = diff([0 this_to_0 T(tt)]);
       good_pre     = isi(1:end-1) > p.min_pre_dur;
       good_post    = isi(2:end) > p.min_post_dur;
       good_t       = this_to_0 > p.min_switch_t & this_to_0 < p.max_switch_t;
       good_switch  = good_pre & good_post & good_t;
        switch_to_0{tt}    = this_to_0(good_switch);

   end
   if ~isempty(this_to_1)
       isi = diff([0 this_to_1 T(tt)]);
       good_pre     = isi(1:end-1) > p.min_pre_dur;
       good_post    = isi(2:end) > p.min_post_dur;
       good_t       = this_to_1 > p.min_switch_t & this_to_1 < p.max_switch_t;
       good_switch  = good_pre & good_post & good_t;
       switch_to_1{tt}    = this_to_1(good_switch);
   end
end

has_switch_to_0 = ~cellfun(@isempty,switch_to_0);
has_switch_to_1 = ~cellfun(@isempty,switch_to_1);
if p.exclude_final
    switch_to_0(has_switch_to_0) = cellfun(@(x) x(1:end-1), switch_to_0(has_switch_to_0),'uniformoutput',0);
    switch_to_1(has_switch_to_1) = cellfun(@(x) x(1:end-1), switch_to_1(has_switch_to_1),'uniformoutput',0);
elseif p.final_only
    switch_to_0(has_switch_to_0) = cellfun(@(x) x(end), ...
        switch_to_0(has_switch_to_0),'uniformoutput',0);
    switch_to_1(has_switch_to_1) = cellfun(@(x) x(end), ...
        switch_to_1(has_switch_to_1),'uniformoutput',0);
end

