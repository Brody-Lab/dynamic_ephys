function [array_data, vec_data, cellid, sessid] = get_behavior_data(datapath, cellid, sessid, p)

disp(['cellid_' num2str(cellid) '_sessid_' num2str(sessid) '_rat_' p.ratname])

try
    % check if this cell has already been analyzed. 
    if ~p.reload & exist([datapath 'cellid_' num2str(cellid) '_sessid_' num2str(sessid) '_rat_' p.ratname '.mat'],'file')==2;
        disp('   loading data...')
        load([datapath 'cellid_' num2str(cellid) '_sessid_' num2str(sessid) '_rat_' p.ratname]);
    else
        % analyze cell
        disp('   computing...')
        [array_data, vec_data, sessid, p.ratname] = package_dyn_phys(cellid);
        [array_data, vec_data]  = cleanup_array_data(array_data, vec_data);
        array_data = compute_state_switches(array_data);
        array_data = compute_gen_state(array_data);
        save([datapath 'cellid_' num2str(cellid) '_sessid_' num2str(sessid) '_rat_' p.ratname], 'array_data','vec_data','sessid','p','cellid')
    end
catch me
    disp(me)
end

