function M = norm_by_peak(M,N)
if nargin < 2
    N = M;
end
    M = M./max(N,[],2);
