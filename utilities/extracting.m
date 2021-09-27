% [d] = extracting(cellbase, colname, selectstr)
% 
% Given a cellbase and a column name(s) and a selctions string, returns a
% cell vector (matrix) containing the indicated column(s), but
% only for the rows of cellbase where selectstr evaluates to true.
%

function [d] = extracting(cellbase, colname, selectstr)

if nargin<3, 
	sel = ones(size(cellbase,1)-1,1);
else
	sel = selecting(cellbase, selectstr);
end;

c = findRnum(cellbase, colname);

d = cellbase(find(sel)+1,c);

