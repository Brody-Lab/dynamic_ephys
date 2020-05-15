function [array_data, vec_data] = cleanup_array_data(array_data, vec_data)

% put all time events relative to stim_start
for i=1:length(array_data)
    array_data(i).left_bups = array_data(i).left_bups - array_data(i).stim_start;
    array_data(i).right_bups = array_data(i).right_bups - array_data(i).stim_start;
    array_data(i).genSwitchTimes = array_data(i).genSwitchTimes - array_data(i).stim_start;   
    array_data(i).spikes = array_data(i).spikes - array_data(i).stim_start;   
    array_data(i).stim_end = array_data(i).stim_end - array_data(i).stim_start;  
    array_data(i).stim_start = array_data(i).stim_start - array_data(i).stim_start;    
end
vec_data.state_0_exits      =  vec_data.state_0_exits   - vec_data.stim_start;
vec_data.state_0_entries    =  vec_data.state_0_entries - vec_data.stim_start;
vec_data.spoke_in           =  vec_data.spoke_in        - vec_data.stim_start;
vec_data.cpoke_out          =  vec_data.cpoke_out       - vec_data.stim_start;
vec_data.cpoke_end          =  vec_data.cpoke_end       - vec_data.stim_start;
vec_data.cpoke_start        =  vec_data.cpoke_start     - vec_data.stim_start;
vec_data.stim_start         =  vec_data.stim_start      - vec_data.stim_start;

