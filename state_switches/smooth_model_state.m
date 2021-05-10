function [array_data] = smooth_model_state(array_data,varargin);
% This function remove several types of possible degenerate model state change predictions
%
% If p.remove_initial_choice, the first state change from no 
% prediction to either L/R is removed
% 
% If p.clear_bad_strengths, then state changes with a strength less 
% than p.bad_strength are removed
% 
% This function does NOT, but SHOULD deal with the following:
% After removing some changes for the above reasions, then we could be left
% with two changes that both go in the same direction. This should probably
% just be replaced with one change at either one of the time points, 
% or perhaps the mean.
%
% If a single click causes a state change, and the state then reverses 
% again, the current method will call that two state changes, when we 
% might want to ignore that. I think it should actually be handled by 
% changing the quantification method. A good alternative would be to fit a 
% line, rather than look at the mean difference of the trajectory. 
%
% Useful debugging trials. Cell 16898, Trials 6, 15, 17, 34, 65, 81
%
%
p = inputParser();
addParameter(p,'remove_initial_choice',1)
addParameter(p,'clear_bad_strengths',1);
addParameter(p,'bad_strength',0);
addParameter(p,'change_bounds',[0 0]);
addParameter(p,'t_buffers',[0 0]);

parse(p,varargin{:});
p = p.Results;

% Remove start of trial state change "initial choice"

for ii=1:length(array_data)
    array_data(ii).ignore_model_switch_to_0 = [];
    array_data(ii).ignore_model_switch_to_1 = [];
    array_data(ii).ignore_0_strength = [];
    array_data(ii).ignore_1_strength = [];
    
    switch_t =  [array_data(ii).model_switch_to_0 ...
        array_data(ii).model_switch_to_1];
    switch_state =  [zeros(size(array_data(ii).model_switch_to_0 ))...
        ones(size(array_data(ii).model_switch_to_1))];
    strengths =  [array_data(ii).model_switch_to_0_strength...
        array_data(ii).model_switch_to_1_strength];
    [switch_t_sorted sort_ind] = sort(switch_t);
    switch_state_sorted = switch_state(sort_ind);
    strengths_sorted    = strengths(sort_ind);
    baddex = zeros(size(strengths));
    
    % if there are no switches, move to the next trial
    if isempty(switch_state_sorted)
        continue;
    end
    
    if p.remove_initial_choice 
        baddex(1) = 1;
    end
    
    baddex(switch_t_sorted < p.t_buffers(1)) = 1;
    max_t = array_data(ii).stim_end - p.t_buffers(2);
    baddex(switch_t_sorted > max_t) = 1;
            
    if p.clear_bad_strengths
        weak_0 = strengths_sorted > -p.bad_strength &...
            switch_state_sorted==0;
        weak_1 = strengths_sorted < p.bad_strength & ...
            switch_state_sorted==1;
        baddex(weak_0 | weak_1) = 1;
    end
    
    % CLEAR BASED ON BOUND CROSSING
    for ss = 1:length(switch_t_sorted);
        if ss == 1
            t0 = 0;
        else
            t0 = switch_t_sorted(ss-1);
        end
        if ss == length(switch_t_sorted)
            tN = array_data(ii).stim_end;
        else
            tN = switch_t_sorted(ss+1);
        end
        pre_T_ind   = array_data(ii).model_T > t0 & ...
            array_data(ii).model_T < switch_t_sorted(ss);
        post_T_ind  = array_data(ii).model_T > switch_t_sorted(ss) & ...
            array_data(ii).model_T < tN;
        pre_mean    = array_data(ii).model_mean(pre_T_ind);
        post_mean    = array_data(ii).model_mean(post_T_ind);
        if switch_state_sorted(ss) == 1
            % currently in state 1, so look for state 0 change bound before cross
            pre_cross   = sum(pre_mean < p.change_bounds(1)) > 0;
            % look for state 1 change bound after cross
            post_cross  = sum(post_mean > p.change_bounds(2)) > 0;
            
        else
            % currently in state 0, so look for state 1 change bound before
            % cross and state 0 change after
            pre_cross   = sum(pre_mean > p.change_bounds(2)) > 0;
            post_cross  = sum(post_mean < p.change_bounds(1)) > 0;
        end
        if ~pre_cross | ~post_cross
            baddex(ss)   = 1;
        end
    end
    
    baddex_0 = find(baddex(switch_state_sorted == 0));
    baddex_1 = find(baddex(switch_state_sorted == 1));

    
    if ~isempty(baddex_0);
        array_data(ii).ignore_0_strength = ...
            [array_data(ii).ignore_0_strength ...
            array_data(ii).model_switch_to_0_strength(baddex_0)];
        array_data(ii).ignore_model_switch_to_0 = [array_data(ii).ignore_model_switch_to_0...
            array_data(ii).model_switch_to_0(baddex_0)];
        array_data(ii).model_switch_to_0_strength(baddex_0) = [];
        array_data(ii).model_switch_to_0(baddex_0) = [];
    end
    if ~isempty(baddex_1);
        array_data(ii).ignore_1_strength = [array_data(ii).ignore_1_strength ...
            array_data(ii).model_switch_to_1_strength(baddex_1)];
        array_data(ii).ignore_model_switch_to_1 = [array_data(ii).ignore_model_switch_to_1...
            array_data(ii).model_switch_to_1(baddex_1)];
        array_data(ii).model_switch_to_1_strength(baddex_1) = [];
        array_data(ii).model_switch_to_1(baddex_1) = [];
    end
end


