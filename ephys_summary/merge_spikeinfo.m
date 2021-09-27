function [allspikes] = merge_spikeinfo(allspikes, spikeinfo)

    allspikes.sessid  = [allspikes.sessid; spikeinfo.sessid];
    allspikes.channel = [allspikes.channel; spikeinfo.channel];
    allspikes.single  = [allspikes.single; spikeinfo.single];
    allspikes.cellid  = [allspikes.cellid; spikeinfo.cellid];
    allspikes.times   = [allspikes.times(:); spikeinfo.times(:)];
    allspikes.ratname = [allspikes.ratname(:); spikeinfo.ratname(:)];
