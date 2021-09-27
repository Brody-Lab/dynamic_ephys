function cmap = colormapRedBlue(nShades,pwr,high)
% return colormap with a red to white to blue gradient
% cmap = colormapRedBlue(nShades,pwr)
% 
% e.g. if showing an image of correlation values, use
%
% image(theCorrelations)
% set(gca,'clim',[-1 1])
% colormap(colormapRedBlue)
%
% to have positive correlations in red, negative in blue
%
%
%
% arguments (optional, provide [] to skip)
%
% nShades - number of colors in the red and blue halves.  size(cmap,1) = nShades*2 + 1
% pwr - used to change gamma 
%


% use default values is none provided
if ~exist('nShades','var') || isempty(nShades)
    nShades = 50;
end

if ~exist('pwr','var') || isempty(pwr)
    pwr = 1;
end

if ~exist('high','var') || isempty(high)
    high=1;
end


% generate gradient of values
grad = linspace(0,high,nShades+1)' .^ pwr;
grad = grad(1:end-1);
gradF = flipud(grad);

%assemble colormap
cmap = [grad grad ones(nShades,1); ...   red gradient
    high.*ones(1,3); ...                       pure white
    ones(nShades,1) gradF gradF];%       blue gradient


