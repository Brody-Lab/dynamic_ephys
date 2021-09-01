function res = compute_rank1_fgta_approx(res, varargin)
p = inputParser;
addParameter(p, 'which_map', 'fgta_resid')
addParameter(p, 'dvlims', [])
addParameter(p, 'tlims', [])
parse(p, varargin{:});
p = p.Results;

if isstruct(res)
    map = res.(p.which_map);
else
    map = res;
end
dvlims  = p.dvlims;
tlims   = p.tlims;
if isempty(dvlims)
    dvlims = res.dv_axis([1 end]);
end
if isempty(tlims)
    tlims = res.t0s([1 end]);
end
goodx = res.dv_axis >= dvlims(1) &  res.dv_axis <= dvlims(2);
goodt = res.t0s >= tlims(1) & res.t0s <= tlims(end);
%goodt(1) = false;
% SVD analysis
% do rank1 approximation
[u,s,v]     = svd(map(goodt,goodx));
s_squared   = diag(s).^2;
s1          = s(1);
u1          = nan(size(res.t0s))';
u1(goodt)   = u(:,1);
v1          = nan(size(res.dv_axis))';
v1(goodx)   = v(:,1);
map_hat     = u1*s1*v1';

rank1_mt    = u1;
rank1_ra    = v1;
% force average firing rate modulation to be positive
if nanmean(rank1_mt) < 0
    rank1_mt = -rank1_mt;
    rank1_ra = -rank1_ra;
end
k           = max(rank1_mt);
alpha       = 1/range(v1);
beta        = s1/alpha;
rank1_mt_n      = rank1_mt*beta;
rank1_ra_n      = rank1_ra*alpha;
rank1_mt_nm     = rank1_mt ./ k;
rank1_ra_nm     = rank1_ra * k * s1;

% see how much variance is explained by different rank approximations
nranks = min(length(s_squared),5);
rank_var = nan(nranks,1);
for ii = 1:nranks
    rank_var(ii) = sum(s_squared(1:ii))./sum(s_squared);
end
res.u       = u;
res.s       = diag(s);
res.v       = v;
res.u1      = u1;
res.s1      = s1;
res.v1      = v1;
res.map_hat = map_hat;
res.alpha   = alpha;
res.beta    = beta;
res.rank1_mt_n  = rank1_mt_n;
res.rank1_ra_n  = rank1_ra_n;
res.rank1_mt_nm = rank1_mt_nm;
res.rank1_ra_nm = rank1_ra_nm;
res.rank_var    = rank_var;

% warning('off','all')
% svdta   = rank1_ra_n';
% svdtan  = svdta - min(svdta);
% svdtan  = svdtan./max(svdtan);
% x       = res.dv_axis;
% [betas,~,~,sigma,mse] = nlinfit(x,svdtan,@dyn_sig,[0, 1/3]);
% warning('on','all')
% svd_betas    = betas;
% svd_delta    = sqrt(diag(sigma)) * tinv(1-0.05/2,sum(~isnan(x))-4);
% svd_sigmas   = svd_delta;
% svd_slope    = betas(2)/4;
% res.svdtan      = svdtan;
% res.svd_betas   = svd_betas;
% res.svd_sigmas  = svd_sigmas;
% res.svd_slope   = svd_slope;