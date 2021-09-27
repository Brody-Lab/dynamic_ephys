function [C] = color_set(nclasses)
 
%cmap = colormap(winter);  
cmap = dyn_cmap.^0.6;

%   cmap = rainbow_colormap;

   if nclasses==1,
      g = round(0.5*size(cmap,1));
      C = cmap(g,:);
      return;
   end;
   
   C = zeros(nclasses, 3);
   
   for i=1:nclasses,
      g = ((i-1)/(nclasses-1))*(size(cmap,1)-1) + 1;
      g = round(g);
      C(i,:) = cmap(g,:);
   end;
