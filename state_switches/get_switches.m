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
addParameter(p,'t_buffers',[0 0]);
addParameter(p,'which_trials',[]);
addParameter(p,'remove_initial_choice', 1);
addParameter(p,'eval_dt', 1e-3);
addParameter(p,'strength_window', .1);
addParameter(p,'model_smooth_wdw', 100);
addParameter(p,'change_bounds', [0 0]);
addParameter(p,'recompute_switches', 0);


parse(p,varargin{:});
p = p.Results;

array_data  = p.array_data;
vec_data    = p.vec_data;

if isempty(array_data) || isempty(vec_data)
    [array_data, vec_data, this_sessid] = package_dyn_phys(cellid);
else
    this_sessid = bdata('select sessid from cells where cellid={S}',cellid);
    if iscell(this_sessid)
        this_sessid = this_sessid{1};
    end
end

% align events to stimulus onset
array_data = cleanup_array_data(array_data, vec_data);
% add generative state switches to array_data
array_data = compute_state_switches(array_data);
array_data = compute_gen_state(array_data);
% decide either to use the model or generative switches
switch p.which_switch
    case 'model'
        % if array_data wasn't passed in with model switches, compute them
        % using the saved model mean and the relevant model parameters
        has_switches = ~isfield(array_data, 'model_switch_to_0') ||  ~isfield(array_data, 'model_switch_to_1');
        if has_switches | p.recompute_switches
            array_data = compute_model_switches(array_data, this_sessid, ...
                'remove_initial_choice', p.remove_initial_choice,...
                'eval_dt', p.eval_dt, ...
                'strength_window', p.strength_window, ...
                'clear_bad_strengths', p.clear_bad_strengths, ...
                'bad_strength', p.bad_strength, 'fit_line', p.fit_line,...
                'model_smooth_wdw',p.model_smooth_wdw,...
                'change_bounds',p.change_bounds,...
                't_buffers',p.t_buffers);
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
   max_t_buff = T(tt)-p.t_buffers(2);
   if ~isempty(this_to_0)
       isi = diff([0 this_to_0 T(tt)]);
       good_pre     = isi(1:end-1) > p.min_pre_dur;
       good_post    = isi(2:end) > p.min_post_dur;
       good_t       = this_to_0 > p.min_switch_t & ...
           this_to_0 < p.max_switch_t &...
           this_to_0 > p.t_buffers(1) & ...
           this_to_0 < max_t_buff;
       % this is a little clumsy, but we are reimposing the t_buffers for
       % generative switches
       good_switch  = good_pre & good_post & good_t;
        switch_to_0{tt}    = this_to_0(good_switch);

   end
   if ~isempty(this_to_1)
       isi = diff([0 this_to_1 T(tt)]);
       good_pre     = isi(1:end-1) > p.min_pre_dur;
       good_post    = isi(2:end) > p.min_post_dur;
       good_t       = this_to_1 > p.min_switch_t & ...
           this_to_1 < p.max_switch_t & ...
           this_to_1 > p.t_buffers(1) & ...
           this_to_1 < max_t_buff;
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

