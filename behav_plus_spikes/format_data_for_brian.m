

%%
cellid = 18181;
celldata    = dyn_cell_packager(cellid);
%%
ratname = celldata.ratname;
sessid = celldata.sessid;
%%
rat_cells = bdata('select cellid from cells where ratname="{S}"',ratname);
%%
cellids = [17786 17787 17788 17791 17793 17795 17797 17798 17799 17802 ...
    17803 17804 18181 18185 18186 18190 18192 18193 18195 18197 ...
    18204 18363 18466 18467];
%%
sessids = [];
for cc = 1:length(cellids)
    sessids(cc) = bdata('select sessid from cells where cellid={S}',cellids(cc));
end
%%
unique_sessids = unique(sessids);


for ss = 1:length(unique_sessids)
    this_sessid = unique_sessids(ss);
    disp(this_sessid)
    sess_cellid = cellids(sessids == this_sessid);
    clear rawdata
    for cc = 1:length(sess_cellid)
        this_cellid = sess_cellid(cc);
        [array_data, vec_data, ~, ~] = package_dyn_phys(this_cellid);
        nt = length(array_data);
        if cc == 1
            rawdata(nt) = struct(); %#ok<SAGROW>
        end
        
        for tt = 1:nt
            % get spikes relative to state 0 and realign to 
            ref = array_data(tt).left_bups(1);
            rawdata(tt).spike_times{cc} = array_data(tt).spikes - ref; %#ok<SAGROW>
            if cc == 1
                rawdata(tt).leftbups = array_data(tt).left_bups - ref;
                rawdata(tt).rightbups = array_data(tt).right_bups - ref;
                rawdata(tt).T = array_data(tt).stim_end - ref;
                rawdata(tt).pokedR = vec_data.pokedR(tt);
            end
            
        end
    end
    save_dir = '/Volumes/jtb3/projects/pbups_dyn/data/brian_phys_format';
    save_name = sprintf('dyn_phys_%s_%i.mat', ratname, this_sessid);
    save(fullfile(save_dir, save_name),  'rawdata', 'sess_cellid');
end
%%
[fh ax] = example_cell_psth('cells',[17793 17791],'pause',1);

