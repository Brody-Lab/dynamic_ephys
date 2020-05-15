function [M_sorted, i_sorted, j_sorted] = sort_by_peak(M,N)
if nargin < 2
    N = M;
end
[pk, pki]  = max(N,[],2);

[j_sorted, i_sorted] = sort(pki); 

M_sorted = M(i_sorted,:);
