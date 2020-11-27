function [gs] = compute_state(T, switchTimes, endState,dt)
% Returns a vector with the generative state at each point in the trial
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
