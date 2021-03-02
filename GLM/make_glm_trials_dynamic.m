function rawData= make_glm_trials_dynamic(Cells,varargin)
p = inputParser;
p.KeepUnmatched=true;
p.addParameter('ref_event','cpoke_in',@(x)validateattributes(x,{'char'},{'nonempty'}));
p.addParameter('samplingFreq',1e3,@(x)validateattributes(x,{'numeric'},{'scalar'}));
p.addParameter('removeViolations',true,@(x)validateattributes(x,{'logical'},{'scalar'}));
p.addParameter('removeStimTrials',true,@(x)validateattributes(x,{'logical'},{'scalar'}));
p.addParameter('nClickBins',1,@(x)validateattributes(x,{'numeric'},{'scalar','positive'}));
p.addParameter('separate_stereo_click',true,@(x)validateattributes(x,{'logical'},{'scalar'}));
p.addParameter('separate_clicks_by','latency',@(x)validatestring(x,{'latency','time'}));
p.addParameter('which_trials',[]);
p.addParameter('use_generative_switches',1);
p.addParameter('use_model_switches',1);
p.addParameter('peh',[]);
p.addParameter('pd',[]);

% time relative to reference event over which you include spikes (make sure the window
% over which you include spiking data in your input structure is the same or smaller)
p.addParameter('spikeWindowS',[-1.25 4],@(x)validateattributes(x,{'numeric'},{'numel',2}));

p.parse(varargin{:});
params = p.Results;

%{
    Things that need to be defined
    --------------------------
    ?ratname
    ?sess_date
    ?mat_file_names
    ?sessid
    ?violations
    ?laser_on
    ?gammas
    ?stim_dur_s
    stim_dur_s % from Cells.Trials. - always shorter than
    stim_dur_s_actual ... not sure what these two are
    ?total_rate
    ?ncells
    ?is_hit
    ?pokedR
    state_ts struct('cpoke_in','cpoke_out','spoke','left_clicks', ...
        'right_clicks','left_click_trial,'right_click_trial')
        % based on Cells.Trials.stateTimes
    spike_ts struct('cpoke_in','cpoke_out','spoke') % based on
        Cells.spike_time_s - like each field is a cell array with length ncells
            each entry is a cell array with length ntrials
            each entry has time stamps of spikes relative to an event like
            and within some window of that event
 
            this seems like kind of
            a crazy way to do this when you could just have the spike times
            and compute their relation as needed
            
%}
sessid = Cells.sessid;
mat_file_name = Cells.mat_file_name;
events = {'cpoke_out', 'cpoke_in','spoke'};

if isempty(params.pd)
    [ratname, sess_date, pd] = bdata(['select ratname, sessiondate, protocol_data ' ...
        'from sessions where sessid={S}'],sessid);
    pd = pd{1};
else
    pd = params.pd;
    ratname = Cells.ratname;
    sess_date = Cells.sess_date;
end

if isempty(params.peh)
    peh         = get_peh(sessid);
else
    peh = params.peh;
end



if isfield(Cells,'c_ts')
    ncells = length(Cells.c_ts);
    ts = Cells.c_ts;
    cellid = Cells.cellid;
else
    cellid      = bdata('select cellid from cells where sessid={S}',sessid);
    ncells      = length(cellid);
end
%%
clear_bad_strengths = 1;
bad_strength    = 0;
fit_line = 1;
exclude_final = 0;
final_only = 0;
which_switch = 'model';
[model_to_0, model_to_1, array_data, vec_data] = ...
    get_switches(cellid, 'which_switch',which_switch,...
    'clear_bad_strengths', clear_bad_strengths, ...
    'bad_strength', bad_strength, 'fit_line', fit_line,...
    'exclude_final', exclude_final, 'final_only', final_only);

switch_to_0 = {array_data.switch_to_0};
switch_to_1 = {array_data.switch_to_1};

switch_tnums = [array_data.trialnum];
%%
violations  = pd.violations;
laser_on    = pd.stims;
ntrials     = length(violations);
bupsdata    = pd.bupsdata(1:ntrials);
gammas      = cellfun(@(x) x.gamma, bupsdata);
nl          = cellfun(@(x) length(x.left), bupsdata);
nr          = cellfun(@(x) length(x.right), bupsdata);
stim_dur_s  = pd.samples;

total_rate  = round(nanmean((nl+nr)./stim_dur_s));

is_hit      = pd.hits;
pokedR      = (pd.sides=='r') == pd.hits;

state_ts    = get_state_ts(peh,pd);

% spikes
for cc = 1:ncells
    if ~isfield(Cells,'c_ts')
        ts = bdata('select ts from spktimes where cellid={S}',cellid(cc));
    end
    for ee = 1:length(events)
        this_event = events{ee};
        for tt = 1:ntrials
            ref_ts = state_ts.(this_event)(tt);
            min_ts = ref_ts + params.spikeWindowS(1);
            max_ts = ref_ts + params.spikeWindowS(2);
            this_ts = ts{1}(ts{1} >= min_ts & ts{1} <= max_ts);
            spike_ts.(this_event){cc}{tt} = this_ts - ref_ts;
        end
    end
end


rawData.rat = ratname;
rawData.param.rat = ratname;
rawData.param.sess_date = sess_date; % 1 session at a time XXXX-XX-XX
%rawData.sess_time = sess_time; % not actually a field
rawData.mat_file_name = mat_file_name; % where we saved the Cells .mat file
rawData.sessid = sessid;
rawData.param.samplingFreq = params.samplingFreq; % Hz
rawData.param.aligned_to = params.ref_event;
rawData.param.ncells = ncells;

rawData.timings = events;
%rawData.timings = events{2};

% if params.removeViolations
%     exclude_trials = violations;
% end
% if params.removeStimTrials
%     exclude_trials = exclude_trials | laser_on;
% end
% % remove side LED and free choice trials by default
% exclude_trials = exclude_trials | abs(gammas) > 90;
exclude_trials = true(ntrials,1);
exclude_trials([array_data.trialnum]) = false;

%% preallocate structure in memory
duration = round(diff(params.spikeWindowS)*params.samplingFreq);
rawData.trial = struct();
rawData.trial(sum(~exclude_trials)).duration = duration; % preallocate
rawData.nTrials = sum(~exclude_trials);
switch params.separate_clicks_by
    case 'time'
        right = linspace(0,params.samplingFreq*max(stim_dur_s),params.nClickBins+1);
    case 'latency'
        binEdges = [0 expinv((1:params.nClickBins-1)/params.nClickBins,params.samplingFreq/total_rate) Inf];
end
good_trials = find(~exclude_trials);



%%
for t = 1:length(array_data)
    original_t = good_trials(t);
    trial_start_time = state_ts.(params.ref_event)(original_t) + params.spikeWindowS(1);
    rawData.trial(t).duration = duration;
    rawData.trial(t).ishit = is_hit(original_t);
    
    for e=1:length(events)
        rawData.trial(t).(events{e}) = state_ts.(events{e})(original_t) - trial_start_time;
        if isnan(rawData.trial(t).(events{e}))
            rawData.trial(t).(events{e}) = [];
        else
            rawData.trial(t).(events{e}) = round(rawData.trial(t).(events{e})*params.samplingFreq);
        end
    end
    rawData.trial(t).gamma = gammas(original_t);
    rawData.trial(t).stim_dur = stim_dur_s(original_t);
    rawData.trial(t).pokedR = pokedR(original_t);
    for c = 1:ncells
        rawData.trial(t).(['sptrain',num2str(c)]) = round( ...
            (spike_ts.(params.ref_event){c}{original_t} - params.spikeWindowS(1) ) * params.samplingFreq);
    end
    %% clicks (stereo click gets its own covariate, and depending on the value of nclickbins
    this_trial_left_clicks = params.samplingFreq*(...
        state_ts.left_clicks(state_ts.left_click_trial==original_t) - trial_start_time);
    this_trial_right_clicks = params.samplingFreq*(...
        state_ts.right_clicks(state_ts.right_click_trial==original_t) - trial_start_time);
    
    if strcmp(params.separate_clicks_by,'latency')
        all_clicks = [this_trial_left_clicks' this_trial_right_clicks'];
        click_idx = [zeros(1,numel(this_trial_left_clicks)),ones(1,numel(this_trial_right_clicks))];
        [all_clicks,sort_idx] = sort(all_clicks);
        click_idx = click_idx(sort_idx);
        itis = [Inf,diff(all_clicks)];
        if itis(2) == 0
            itis(2)=Inf;
        else
            warning('No stereo click at start of this trial!');
        end
        this_trial_left_latencies = itis(click_idx==0);
        this_trial_right_latencies = itis(click_idx==1);
    end
    if params.separate_stereo_click
        if this_trial_left_clicks(1)==this_trial_right_clicks(1)
            rawData.trial(t).stereo_click = this_trial_left_clicks(1);
            this_trial_left_clicks(1)=[];
            this_trial_right_clicks(1)=[];
            if strcmp(params.separate_clicks_by,'latency')
                this_trial_left_latencies(1)=[];
                this_trial_right_latencies(1)=[];
            end
        else
            warning('No stereo click at start of this trial!');
        end
        if t==1
            rawData.timings = union(rawData.timings,{'stereo_click'});
        end
    end
    for i=1:params.nClickBins
        if params.nClickBins==1
            idx_str = '';
        else
            idx_str = num2str(i);
        end
        if i==1
            if strcmp(params.separate_clicks_by,'time')
                if params.separate_stereo_click
                    these_bin_edges = binEdges + rawData.trial(t).stereo_click;
                else
                    these_bin_edges = binEdges + min([this_trial_left_clicks(:);this_trial_right_clicks(:)]);
                end
            end
        end
        
        if t==1
            rawData.timings = union(rawData.timings,{['left_clicks',idx_str] ['right_clicks',idx_str] ['all_clicks',idx_str]});
        end
        switch params.separate_clicks_by
            case 'time'
                rawData.trial(t).(['left_clicks',idx_str])  = ...
                    this_trial_left_clicks(this_trial_left_clicks>=these_bin_edges(i) ...
                    & this_trial_left_clicks<these_bin_edges(i+1)) ;
                rawData.trial(t).(['right_clicks',idx_str])  = ...
                    this_trial_right_clicks(this_trial_right_clicks>=these_bin_edges(i) ...
                    & this_trial_right_clicks<these_bin_edges(i+1)) ;
                rawData.trial(t).(['all_clicks',idx_str])  = ...
                    union(rawData.trial(t).(['left_clicks',idx_str]),rawData.trial(t).(['right_clicks',idx_str])) ;
            case 'latency'
                rawData.trial(t).(['left_clicks',idx_str])  = ...
                    this_trial_left_clicks(this_trial_left_latencies>=binEdges(i) ...
                    & this_trial_left_latencies<binEdges(i+1)) ;
                rawData.trial(t).(['right_clicks',idx_str])  = ...
                    this_trial_right_clicks(this_trial_right_latencies>=binEdges(i) ...
                    & this_trial_right_latencies<binEdges(i+1)) ;
                rawData.trial(t).(['all_clicks',idx_str])  = ...
                    union(rawData.trial(t).(['left_clicks',idx_str]),rawData.trial(t).(['right_clicks',idx_str])) ;
        end
        
        
    end
    if params.use_generative_switches
        % add timings for switches
        if t==1
            rawData.timings = union(rawData.timings,{['switch_to_0'] ['switch_to_1'] }); %#ok<NBRAK>
        end
        
        switch_t = t;%find(switch_tnums==original_t);
        if ~isempty(switch_t) 
            off = +rawData.trial(t).stereo_click*1e-3 - array_data(t).right_bups(1);

            this_to_0 = sort(switch_to_0{switch_t}) + ...
                off;
            this_to_1 = sort(switch_to_1{switch_t}) + ...
                off;
            
            if isempty(this_to_0)
                rawData.trial(t).switch_to_0 = [];
            else
                rawData.trial(t).switch_to_0 = round(params.samplingFreq*(this_to_0));
            end
            if isempty(this_to_1)
                rawData.trial(t).switch_to_1 = [];
            else
                rawData.trial(t).switch_to_1 = round(params.samplingFreq*(this_to_1));
            end
        end
    end
    if params.use_model_switches
        % add timings for switches
        if t==1
            rawData.timings = union(rawData.timings,{['model_to_0'] ['model_to_1'] }); %#ok<NBRAK>
        end
        
        switch_t = t;%find(switch_tnums==original_t);
        if ~isempty(switch_t)
            off = +rawData.trial(t).stereo_click*1e-3 - array_data(t).right_bups(1);
            this_to_0 = sort(model_to_0{switch_t}) + ...
                off;
            this_to_1 = sort(model_to_1{switch_t}) + ...
                off ;
            
            if isempty(this_to_0)
                rawData.trial(t).model_to_0 = [];
            else
                rawData.trial(t).model_to_0 = round(params.samplingFreq*(this_to_0));
            end
            if isempty(this_to_1)
                rawData.trial(t).model_to_1 = [];
            else
                rawData.trial(t).model_to_1 = round(params.samplingFreq*(this_to_1));
            end
        end
    end
    
end

end
%%
