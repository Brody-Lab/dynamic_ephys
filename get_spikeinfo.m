function spikeinfo = get_spikeinfo(ratname)

% get session ids for this rat and date range
sessid   = bdata('select sessid from sessions where ratname="{S}"',ratname);
    
% remove sessions without phys
which_phys = false(size(sessid));
for ss = 1:length(sessid)
    this_sessid = sessid(ss);
    is_phys_sess = ~isempty(bdata('select sessid from phys_sess where ratname="{S}" and sessid={S}',ratname,num2str(this_sessid)));
    if is_phys_sess
        which_phys(ss) = true;
    end
end
phys_sessid = sessid(which_phys);
%%

% select the cellids from the phys sessions
cellid      = [];
is_single   = [];
channels    = [];
sessid      = [];
for i=1:length(phys_sessid)
[t_cellid, t_is_single, t_channels,t_sessid] = bdata(['select cellid, single, sc_num,sessid from cells where sessid="{S}"'], num2str(phys_sessid(i)));
cellid = [cellid; t_cellid];
is_single = [is_single; t_is_single];
channels = [channels; t_channels];
sessid = [sessid; t_sessid];
end

%%
% get the spike times for these cells
spktime = cell(1,length(cellid));
for cc = 1:length(cellid)
    this_ts = bdata('select ts from spktimes where cellid={S}',num2str(cellid(cc)));
    spktime{cc} = this_ts{:};
    
end

spikeinfo.times = spktime;
spikeinfo.cellid = cellid;
spikeinfo.single = is_single;
spikeinfo.channel = channels;
spikeinfo.sessid  = sessid;

