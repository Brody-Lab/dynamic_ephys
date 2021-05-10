function array_data = get_model_state_durs(array_data);
for ii = 1:length(array_data)
    m_switch_t =  sort([array_data(ii).model_switch_to_0 ...
        array_data(ii).model_switch_to_1]);
    
    tN = array_data(ii).stim_end;
    array_data(ii).model_state_durs = diff([0 m_switch_t tN]);
    array_data(ii).model_switches_per_time  = length(m_switch_t) / tN;
    
    array_data(ii).gen_state_durs = ...
        diff([0 array_data(ii).genSwitchTimes tN]);
    array_data(ii).gen_switches_per_time = ...
        length(array_data(ii).genSwitchTimes) / tN;
end