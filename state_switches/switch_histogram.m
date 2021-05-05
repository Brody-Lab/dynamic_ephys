cellid = 18181;
which_switch = 'generative'
[switch_to_0, switch_to_1, array_data, vec_data] = ...
    get_switches(cellid, 'which_switch',which_switch,...
    'clear_bad_strengths', 0, ...
    'exclude_final', 0, 'final_only', 0);
figure; histogram(cellfun(@length,switch_to_1))

%%
both_switches = [switch_to_0{:} switch_to_1{:}]
T = {array_data.stim_end}
diffs = cellfun(@(x,y) x-y, T, switch_to_0,...
    'uniformoutput',0,'errorfun',@(x) nan);
histogram(diffs,[0:.25:2])