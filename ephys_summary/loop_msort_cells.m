fpath = '/Volumes//brody/physdata/Ahmed_SpikeGadgets/';
assert(exist(fpath)>0)
fn = 'spikes_@PBupsDyn_Ahmed_H191_190624a.mat';
fprintf('loading %s',fn)
tic
d= load(fullfile(fpath,fn))
toc
%%
fofwv = [d.wave{:}];
[~, ordering] = sort(max(fofwv))
 inc = 1;
for ww = 1:length(d.wave);
    figure(1); clf;
    for pp = 1:25;
        inc = inc +1;
        subplot(5,5,pp);
        plot(d.wave{ordering(inc)}')
    end;
    pause;
end