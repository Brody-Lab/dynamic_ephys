function cmap = colormapLinear(max_color,nShades,min_color,pwr)
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

if ~exist('min_color','var') || isempty(min_color)
    min_color = [1 1 1];
end

% generate gradient of values
grad1 = linspace(max_color(1),min_color(1),nShades+1)' .^ pwr;
grad2 = linspace(max_color(2),min_color(2),nShades+1)' .^ pwr;
grad3 = linspace(max_color(3),min_color(3),nShades+1)' .^ pwr;


% assemble colormap
cmap = flipud([grad1 grad2 grad3]);

