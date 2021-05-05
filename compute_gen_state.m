function [array_data] = compute_gen_state(array_data)
% adds a field to array_data "gen_state" which is a vector of 0s or 1s to indicate state
% adds a field to array_data "gen_state_duration" which is a vector of how long since the last state change.
% all state changes are in terms of the generative state
% all state changes *do not* take silent transitions into account


for i=1:length(array_data)
    T = (array_data(i).stim_end - array_data(i).stim_start);
    ns = round((array_data(i).stim_end - array_data(i).stim_start)/1e-4);
    array_data(i).gen_state = compute_gen_state_trial(T, array_data(i).genSwitchTimes-array_data(i).stim_start, array_data(i).genEndState, 1e-4);
   
    
    trial_t = 1e-4*ones(size(array_data(i).gen_state(1:end)));
    gen_state_dur       = cumsum(trial_t);
    gen_state_switches  = [0; diff(array_data(i).gen_state)];
    gen_state_switch_ix = find(gen_state_switches);
    gen_state_diff      = zeros(size(array_data(i).gen_state));
    
    for ss = 1:length(gen_state_switch_ix)
        ii = gen_state_switch_ix(ss);
        gen_state_dur(ii:end) = gen_state_dur(ii:end) - gen_state_dur(ii);
    end
    array_data(i).gen_state_duration = gen_state_dur;
end





