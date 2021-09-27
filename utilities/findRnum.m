%findRnum   [n] = findRnum(R, str)  Find number of given desc(s).
%
% Returns the number (that is, the order indexing in an R-file) of
% a particular descriptor. For example, the typical index for
% 'class' might be 1, indicating that in a trial cell, the first
% cell contains the class number. Returns empty when it can't find 
% the marker. 
%
% If str is a cell array of strings, does this for all the
% elements. 
%
% Class of returned n is double. Will have same size as the str
% cell array.
%


function [n] = findRnum(R, str)
   
   if ~iscell(str),
      n = find(strcmp(str, R(1,:)));

	  if isempty(n),
		  fprintf(2, ['Warning: couldn''t find marker ' str ...
			  ' in descriptor head\n']);
	  end;
	  return;
   end;
   
   
   n = zeros(size(str));
   
   for i=1:rows(str),
	   for j=1:cols(str)
		   gu = find(strcmp(str{i,j}, R(1,:)));
		   
		   if isempty(gu),
			   fprintf(2, ['Warning: couldn''t find marker ' ...
				   str{i,j} ' in descriptor\n']);
			   n(i, j) = 0;
		   else
			   n(i,j) = gu;
		   end;
	   end;
   end;
