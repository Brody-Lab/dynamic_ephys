function cmap = dyn_cmap(nShades)
lo = 0;
hi = .549;
mid = .2118;
mygreen = [lo  hi  mid];
myblue =  [lo  mid hi];


if ~exist('nShades','var') || isempty(nShades)
    nShades = 50;
end


% generate gradient of values
lowToHi = linspace(lo,hi,nShades+1)' ;
lowToMid = linspace(lo,mid,nShades+1)' ;
midToHi = linspace(mid,hi,nShades+1)' ;



%assemble colormap
cmap = [lowToMid midToHi hi*ones(nShades+1,1); ...   blue gradient
    flipud(lowToMid) hi*ones(nShades+1,1)  flipud(midToHi)];%       blue gradient
