function [array_data] = compute_accumulation_state_switches(array_data);

for i=1:length(array_data)
    left = array_data(i).left_bups;
    right = array_data(i).right_bups;
    clicks = [right' left'; ones(size(right))' -ones(size(left))'];
    [~,dex] = sort(clicks(1,:),2);
    clicks = clicks(:,dex);
    clickdiff = cumsum(clicks(2,:));
    switches_to_1 = find(clickdiff(3:end) == 0 & clickdiff(2:end-1) < 0)+2;
    switches_to_0 = find(clickdiff(3:end) == 0 & clickdiff(2:end-1) > 0)+2;


    array_data(i).accumulation_switch_to_0 = clicks(1,switches_to_0); 
    array_data(i).accumulation_switch_to_1 = clicks(1,switches_to_1); 
    if isempty(switches_to_1);
    array_data(i).accumulation_switch_to_1 = []; 
    end
    if isempty(switches_to_0);
    array_data(i).accumulation_switch_to_0 = []; 
    end
end

