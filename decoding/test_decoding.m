cellid=bdata('select cellid from cells where sessid=565672');
%%
%%
cellid = [17791       17793       17795       17797       17798       17799       17802       17803       17804];
example_cell_psth('cells',cellid)
%%
cellid = [18185       18186];
example_cell_psth('cells',cellid)
%%
sessid  = 599049;
cellid = [17791       17793       17795       17797       17798       17799       17802       17803       17804];
cellid = [18185       18186];
cellid = [17786       17787       17788];
cellid = bdata('select cellid from cells where sessid=511311')
[array_data, vec_data] = package_dyn_phys(cellid(1));


ncells  = length(cellid);

ntrials = length(array_data);

tx = [-1:.1:0];
y = vec_data.pokedR;

ntimes = length(tx) - 1;

n = min(sum(y), sum(~y));
ntrials = 2*n;
left_ind    = find(y==0);
right_ind   = find(y==1);
trialset    = sort([left_ind(randperm(n)); right_ind(randperm(n))]);
X = nan(ntrials, ncells,ntimes);
hits = nan(ntrials,ntimes);

for cc = 1:ncells
    [array_data, vec_data] = package_dyn_phys(cellid(cc));
    
    for ti = 1:ntrials
        for ii = 1:ntimes
            tt= trialset(ti);
            
            this_spikes = array_data(tt).spikes;
            
            t0 = array_data(tt).stim_end+tx(ii);
            tn = array_data(tt).stim_end+tx(ii+1);
            this_spk_ct = sum(this_spikes > t0 & this_spikes < tn);
            
            X(ti,cc,ii) = this_spk_ct;
            
        end
    end
end

for ii = 1:ntimes
    thisX = X(:,:,ii);
    [b,info] = lassoglm(thisX, y(trialset),  'binomial');
    t = 1;
    b = [ info.Intercept(t); b(:,t)];
    yy = glmval(b, thisX, 'logit');
    yhat = yy > .5;
    hits(:,ii) = yhat == y(trialset);
end

figure; plot(tx(1:end-1),mean(hits))
