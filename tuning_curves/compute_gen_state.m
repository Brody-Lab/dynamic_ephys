function [array_data] = compute_gen_state(array_data)
% adds a field to array_data "gen_state" which is a vector of 0s or 1s to indicate state
% adds a field to array_data "gen_state_duration" which is a vector of how long since the last state change.
% all state changes are in terms of the generative state
% all state changes *do not* take silent transitions into account


for i=1:length(array_data)
    T = (array_data(i).stim_end - array_data(i).stim_start);
    ns = round((array_data(i).stim_end - array_data(i).stim_start)/1e-4);
    array_data(i).gen_state = compute_gen_state_trial(T, array_data(i).genSwitchTimes-array_data(i).stim_start, array_data(i).genEndState, 1e-4);
   
    count = 1e-4;
    array_data(i).gen_state_duration(1) = 1e-4;
    for j=2:length(array_data(i).gen_state)
        if array_data(i).gen_state(j) == array_data(i).gen_state(j-1)
            count = count + 1e-4;
        else
            count = 1e-4;
        end
        array_data(i).gen_state_duration(j) = count;
    end
end





