function [bindex] = bin_aval(avals, this_val)
    centers     = avals;
    edges       = [-Inf; centers(1:end-1) + diff(centers)/2; Inf];
    [n,edges]   = histcounts(this_val,edges);
    bindex      = find(n);
end
