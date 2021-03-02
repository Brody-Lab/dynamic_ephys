function res = compute_rank1_decomp(map)

% SVD analysis
% do rank1 approximation
[u,s,v]     = svd(map);
s_squared   = diag(s).^2;
s1          = s(1);
u1          = u(:,1);
v1          = v(:,1);
map_hat     = u1*s1*v1';
alpha       = 1/range(v1);
beta        = s1/alpha;
rank1_mt    = u1;
rank1_ra    = v1;
% force average firing rate modulation to be positive
if mean(rank1_mt) < 0
    rank1_mt = -rank1_mt;
    rank1_ra = -rank1_ra;
end
k = max(rank1_mt);
rank1_mt_n      = rank1_mt*beta;
rank1_ra_n      = rank1_ra*alpha;
rank1_mt_nm     = rank1_mt ./ k;
rank1_ra_nm     = rank1_ra * k * s1;


res.rank1_mt    = rank1_mt;
res.rank1_ra    = rank1_ra;
res.rank1_mt_n  = rank1_mt_n;
res.rank1_ra_n  = rank1_ra_n;
res.rank1_mt_nm = rank1_mt_nm;
res.rank1_ra_nm = rank1_ra_nm;