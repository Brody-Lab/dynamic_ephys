function [array_data] = compute_state_switches(array_data)
% computes a list of switch to state 0, and a list of switch to state 1
% all state transitions are generative state
% does not account for silent transitions

for i=1:length(array_data)
    t = array_data(i).genSwitchTimes;
    switch_to_0 = [];
    switch_to_1 = [];

    if array_data(i).genEndState
        dex = 1;
    else
        dex = 0;
    end
    for j = length(t):-1:1
        if dex
            switch_to_1 = [switch_to_1 t(j)];
        else
            switch_to_0 = [switch_to_0 t(j)];       
        end
        dex = ~dex;
    end
    array_data(i).switch_to_0 = switch_to_0;
    array_data(i).switch_to_1 = switch_to_1;
    
end

