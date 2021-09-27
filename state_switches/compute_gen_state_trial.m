function [gs] = compute_gen_state_trial(T, switchTimes, endState,dt)
% For a single trial for duration T, with switchTimes, and an endState, computes the state value of each point in time.
%
% Returns a vector with the generative state at each point in the trial
%
% Tells you the generative state, does not account for silent transitions

    switchTimes = switchTimes(switchTimes >=dt);
    switches = zeros(length(dt:dt:T),1);
    switches(round(switchTimes*(1/dt))) = 1;
    culm = mod(cumsum(switches),2); 
    if culm(end) == endState
        gs = culm;
    else
        gs = ~culm;
    end
end

