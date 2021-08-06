function [array_data] = compute_model_state(array_data,model_mean,change_bound)
%builds fields, model_state, and model_state_duration, model_switch_to_0/1

    % model_state = sign(model trajectory)
    % model_switch_to_0 = 
    % model_switch_to_1 = 
    % model_mean = model mean trajectory
if nargin  < 3
    change_bound = 0;
end

if length(array_data) ~= length(model_mean)
    error('array_data and model_mean are different sizes') 
end

for i=1:length(array_data)
    temp = sort([array_data(i).left_bups; array_data(i).right_bups]);
    if length(temp) >= 3
        firstclk = temp(3);
        clkdex = ceil(firstclk/1e-3);
    else
        firstclk = array_data(i).left_bups(1);
        clkdex = ceil(firstclk/1e-3);
    end
    array_data(i).model_mean = model_mean(i).mean;
    array_data(i).model_mean(1:clkdex) = 0;
    array_data(i).model_state = double(array_data(i).model_mean > change_bound);
    array_data(i).model_state(1:clkdex) = 0.5;
    array_data(i).model_T = model_mean(i).T;
    array_data(i).model_switch_to_0 = model_mean(i).T(find(diff(array_data(i).model_state) < 0));
    array_data(i).model_switch_to_1 = model_mean(i).T(find(diff(array_data(i).model_state) > 0));

end





