function [array_data] = smooth_model_state(array_data,p);
% This function remove several types of possible degenerate model state change predictions
%
% If p.remove_initial_choice, the first state change from no prediction to either L/R is removed
% 
% If p.clear_bad_strengths, then state changes with a strength less than p.bad_strength are removed
% 
% This function does NOT, but SHOULD deal with the following:
% After removing some changes for the above reasions, then we could be left with two changes that both go in the same direction. This should probably just be replaced with one change at either one of the time points, or perhaps the mean.
%
% If a single click causes a state change, and the state then reverses again, the current method will call that two state changes, when we might want to ignore that. I think it should actually be handled by changing the quantification method. A good alternative would be to fit a line, rather than look at the mean difference of the trajectory. 
%
% Useful debugging trials. Cell 16898, Trials 6, 15, 17, 34, 65, 81
%
%

% Remove start of trial state change "initial choice"
if p.remove_initial_choice
for i=1:length(array_data)
    if (length(array_data(i).model_switch_to_0) >= 1) & (length(array_data(i).model_switch_to_1) >= 1)
        % had multiple state changes, need to determine which one was first
        if array_data(i).model_switch_to_0(1) < array_data(i).model_switch_to_1(1)
            array_data(i).ignore_model_switch_to_0 = array_data(i).model_switch_to_0(1);
            array_data(i).model_switch_to_0(1) = [];     
            array_data(i).ignore_0_strength = array_data(i).model_switch_to_0_strength(1);
            array_data(i).model_switch_to_0_strength(1) = [];      
        else
            array_data(i).ignore_model_switch_to_1 = array_data(i).model_switch_to_1(1);
            array_data(i).model_switch_to_1(1) = [];           
            array_data(i).ignore_1_strength = array_data(i).model_switch_to_1_strength(1);
            array_data(i).model_switch_to_1_strength(1) = []; 
        end
    elseif length(array_data(i).model_switch_to_0) == 0
        % Only had one initial choice, 
        array_data(i).ignore_model_switch_to_1 = array_data(i).model_switch_to_1;
        array_data(i).model_switch_to_1 = [];
        array_data(i).ignore_1_strength = array_data(i).model_switch_to_1_strength;
        array_data(i).model_switch_to_1_strength = []; 
    else
        % Only had one initial choice, 
        array_data(i).ignore_model_switch_to_0 = array_data(i).model_switch_to_0;
        array_data(i).model_switch_to_0 = [];
        array_data(i).ignore_0_strength = array_data(i).model_switch_to_0_strength;
        array_data(i).model_switch_to_0_strength = [];      
    end

end
end

% Remove weak changes
if p.clear_bad_strengths 
for i=1:length(array_data)
    baddex = [];
    for j=1:length(array_data(i).model_switch_to_0);
        if array_data(i).model_switch_to_0_strength(j) > -p.bad_strength
            baddex = [baddex j];
        end
    end
    if ~isempty(baddex);
    array_data(i).ignore_0_strength = [array_data(i).ignore_0_strength array_data(i).model_switch_to_0_strength(baddex)];
    array_data(i).ignore_model_switch_to_0 = [array_data(i).ignore_model_switch_to_0 array_data(i).model_switch_to_0(baddex)];
    array_data(i).model_switch_to_0_strength(baddex) = [];
    array_data(i).model_switch_to_0(baddex) = [];
    end

    baddex = [];
    for j=1:length(array_data(i).model_switch_to_1);
        if array_data(i).model_switch_to_1_strength(j) < p.bad_strength
            baddex = [baddex j];
        end
    end
    if ~isempty(baddex)
    array_data(i).ignore_1_strength = [array_data(i).ignore_1_strength array_data(i).model_switch_to_1_strength(baddex)];
    array_data(i).ignore_model_switch_to_1 = [array_data(i).ignore_model_switch_to_1 array_data(i).model_switch_to_1(baddex)];
    array_data(i).model_switch_to_1_strength(baddex) = [];
    array_data(i).model_switch_to_1(baddex) = [];
    end
end
end

% Remove sequential matching state changes?
