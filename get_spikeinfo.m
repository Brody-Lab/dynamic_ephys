function spikeinfo = get_spikeinfo(ratname)

% get list of sessions and cells for this rat
[sessid, cellid] = bdata('select sessid, cellid from cells where ratname="{S}"',ratname);

% get the spike times for these cells
is_single   = false(length(cellid),1);
channels    = zeros(length(cellid),1);
spktime     = cell(1,length(cellid));
for cc = 1:length(cellid)
    this_ts = bdata('select ts from spktimes where cellid={S}',cellid(cc));
    spktime{cc} = this_ts{:};
    [is_single(cc), channels(cc)] = bdata(['select single, sc_num from cells where cellid={S}'], cellid(cc));
end

spikeinfo.times = spktime;
spikeinfo.cellid = cellid;
spikeinfo.single = is_single;
spikeinfo.channel = channels;
spikeinfo.sessid  = sessid;

