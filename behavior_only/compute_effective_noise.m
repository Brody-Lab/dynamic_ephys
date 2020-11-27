function [fits] = compute_effective_noise(fits,data)

% Set up parameters
if length(fits.final) == 9
    phi = fits.final(6);
    tau = fits.final(7);
    s2  = fits.final(3);
elseif length(fits.final) == 7
    phi = fits.final(5);
    tau = fits.final(6);
    s2  = fits.final(3);
elseif length(fits.final) == 8
    phi = fits.final(5);
    tau = fits.final(6);
    s2  = fits.final(3);
elseif length(fits.final) == 6
    phi = fits.final(3);
    tau = fits.final(4);
    s2  = fits.final(2);
end

% compute effective click weight for each click
w  = [];
nt = [];
for i=1:length(data);
    [cl, cr]= make_adapted_cat_clicks(data(i).leftbups, data(i).rightbups, phi, tau);
    w       = [w  mean([cl cr])];
    nt      = [nt length([cl cr])];
end
ab = sum(w.*nt)/sum(nt);

% converting to per click basis, only if using the output of Bing's model
if length(fits.final) == 9
    s2    = s2/40;
end

% Compute final noise
fits.n = .5.*(1+erf(-1*ab./sqrt(s2.*ab.*2)));

numtrials = length(w);

% Bootstrap for confidence intervals
nboots = 100;
n =[];
for i=1:nboots
    disp(i)
    if length(fits.final) == 9
        bphi = randn(1)*fits.se(6) + fits.final(6);
        btau = randn(1)*fits.se(7) + fits.final(7);
        bs2  = randn(1)*fits.se(3) + fits.final(3);
    elseif length(fits.final) == 7
        bphi = randn(1)*fits.se(5) + fits.final(5);
        btau = randn(1)*fits.se(6) + fits.final(6);
        bs2  = randn(1)*fits.se(3) + fits.final(3);
    elseif length(fits.final) == 8
        bphi = randn(1)*fits.se(5) + fits.final(5);
        btau = randn(1)*fits.se(6) + fits.final(6);
        bs2  = randn(1)*fits.se(3) + fits.final(3);
    elseif length(fits.final) == 6
        bphi = randn(1)*fits.se(3) + fits.final(3);
        btau = randn(1)*fits.se(4) + fits.final(4);
        bs2  = randn(1)*fits.se(2) + fits.final(2);
    end
    bs2(bs2 <0)    = 0;
    bphi(bphi < 0) = 0;
    btau(btau < 0) = 0;

%    w  = [];
%    nt = [];
    w = zeros(1,numtrials);
    nt = zeros(1,numtrials);
    for i=1:length(data);
        [cl, cr]= make_adapted_cat_clicks(data(i).leftbups, data(i).rightbups, bphi, btau);

%        w       = [w  mean([cl cr])];
%        nt      = [nt length([cl cr])];
        w(i) = mean([cl cr]);
        nt(i) = length([cl cr]);
    end
    ab = sum(w.*nt)/sum(nt);
    if length(fits.final) == 9
        bs2    = bs2/40;
    end
    n = [n 0.5.*(1+erf(-1*ab./sqrt(bs2.*ab.*2)))];;
end
fits.n_std = std(n);
fits.n_boot = mean(n);

