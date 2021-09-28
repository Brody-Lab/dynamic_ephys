
close all;
clear all;
dp = set_dyn_path(1);

%%
% parameters of simulation
nt = 10;
nx = 21;
ns = 2000000;
END_SIGN_POSITIVE = true;
DEMONSTRATE_WEIGHTING_BAD = false;

% Set up storage
xvec = -floor(nx/2):floor(nx/2);
if END_SIGN_POSITIVE
    END_DEX = xvec > 0;
else
    END_DEX = xvec < 0;
end
nd = sum(END_DEX); % how many delta-function solutions we have to compute
tvec = 1:nt;
s  = zeros(ns,nt);
sf = zeros(nx,nt);
sb = zeros(nx,nt);
f  = zeros(nx,nt);
b  = zeros(nx,nt);
bd = zeros(nx,nt,nd);

% Initialize forward and backwards pass
f(floor(nx/2)+1,1)  = 1; % initialize at 0
b(END_DEX,nt)       = 1; 
b(:,nt)             = b(:,nt)./sum(b(:,nt));
bd(END_DEX,nt,:)    = eye(nd);

% make transition matrix
T    = zeros(nx,nx); 
temp = diag((1/3).*ones(nx,1),-1) + diag((1/3).*ones(nx,1),+1); 
T    = temp(1:nx,1:nx); 
T    = T+eye(nx).*(1/3);

% compute particle simulation forwards
step  = rand(ns,nt);
fstep = step < (1/3);
bstep = (step > (1/3)) & (step <= (2/3));
nstep = ~bstep & ~fstep;
for i=1:nt-1
    s(fstep(:,i),i+1) = s(fstep(:,i),i) + 1;
    s(bstep(:,i),i+1) = s(bstep(:,i),i) - 1;
    s(nstep(:,i),i+1) = s(nstep(:,i),i);
end

% convert particle simulation into distributions
for i=1:nt
    counts = hist(s(:,i),xvec);
    sf(:,i) = counts;
end
sf = sf ./ ns;

% compute backwards particle distribution
if END_SIGN_POSITIVE
    dex = s(:,end) > 0;
else
    dex = s(:,end) < 0;
end
for i=1:nt
    counts = hist(s(dex,i),xvec);
    sb(:,i)= counts;
end
sb = sb./sum(dex);

% Compute forward pass
for i=1:nt-1
    f(:,i+1) = T*f(:,i);
end

% Compute Independent backwards propagation
for i=1:nt-1
    b(:,end-i) = T*b(:,end-i+1);
    for j=1:nd
        bd(:,end-i,j) = T*squeeze(bd(:,end-i+1,j));
    end
end

% Compute Backwards Pass
R = f.*b;
if DEMONSTRATE_WEIGHTING_BAD
    temp = find(END_DEX); %% BAD JUST FOR EXPLORE
end
if DEMONSTRATE_WEIGHTING_BAD
    Rd = f.*squeeze(bd(:,:,1)).*f(temp(1),end); %% BAD JUST FOR EXPLORE
else
    Rd = f.*squeeze(bd(:,:,1));
end
for i=2:nd
    if DEMONSTRATE_WEIGHTING_BAD
        Rd = Rd + f.*squeeze(bd(:,:,i)).*f(temp(i),end); %% BAD JUST FOR EXPLORE
    else
        Rd = Rd + f.*squeeze(bd(:,:,i));
    end
end

% normalize each timestep
for i=1:nt
    R(:,i)  = R(:,i)./sum(R(:,i));
    Rd(:,i) = Rd(:,i)./sum(Rd(:,i));
end
%%
% plot
cax = [0 .4];
fh = figure(1); clf;
fht = 2.5 * 2.5;
fw  = 2.5 * 3.75;
set(fh, 'position', [5 10 fw fht])
plot_old_panels = 0;
if plot_old_panels
    nrow = 3;
else
    nrow = 2;
end

subplot(nrow,3,1);
imagesc(tvec,xvec,f)
title('Forward Distribution, a_0=0')
ylabel('a value')
xlabel('time')
caxis(cax)

subplot(nrow,3,2);
imagesc(tvec,xvec,squeeze(bd(:,:,2)))
title('Backward Distribution, a_N=2')
ylabel('a value')
xlabel('time')
caxis(cax)

subplot(nrow,3,3);
imagesc(tvec,xvec,sf)
title('Forward Particle')
ylabel('a value')
xlabel('time')
caxis(cax)

if plot_old_panels
    subplot(nrow,3,4);
    imagesc(tvec,xvec,b)
    title('Backward Distribution')
    ylabel('a value')
    xlabel('time')
    caxis(cax)
    
    subplot(nrow,3,5);
    imagesc(tvec,xvec,squeeze(bd(:,:,7)))
    title('Posterior Distribution, a=7')
    ylabel('a value')
    xlabel('time')
    caxis(cax)

    subplot(nrow,3,6);
    hold on;
    plot(xvec,sb(:,5),'bo-')
    plot(xvec,Rd(:,5),'kx-')
    plot(xvec,R(:,5),'r--')
    title('Slice at t=5')
    ylabel('Probability')
    xlabel('a value')
    
    subplot(nrow,3,7);
    imagesc(tvec,xvec,R)
    title('Posterior - Analytical')
    ylabel('a value')
    xlabel('time')
    caxis(cax)
    
    subplot(nrow,3,8);
    imagesc(tvec,xvec,Rd)
    title('Posterior - Delta')
    ylabel('a value')
    xlabel('time')
    caxis(cax)
    
    subplot(nrow,3,9);
    imagesc(tvec,xvec,sb)
    title('Posterior Particle')
    ylabel('a value')
    xlabel('time')
    caxis(cax)
else
    subplot(nrow,3,4);
    hold on;
    plot(xvec,sb(:,5),'bo','markersize',8)
    plot(xvec,Rd(:,5),'kx','markersize',5)
    plot(xvec,R(:,5),'-','color',dp.model_color.^2)
    title('Slice at t=5')
    ylabel('Probability')
    xlabel('a value')
    
    subplot(nrow,3,5);
    imagesc(tvec,xvec,R)
    title('Posterior Distribution')
    ylabel('a value')
    xlabel('time')
    caxis(cax)
    
    subplot(nrow,3,6);
    imagesc(tvec,xvec,sb)
    title('Posterior Particle')
    ylabel('a value')
    xlabel('time')
    caxis(cax)
end
end_color = dp.model_color;
colormap(colormapLinear(end_color).^2)

set(findall(fh,'type','axes'),'Box','off');

print(fh,fullfile(dp.fig_dir,'particle_dist_supp'),'-dsvg','-painters')

if DEMONSTRATE_WEIGHTING_BAD
    disp('WARNING, THESE RESULTS SHOULD BE WRONG')
end



