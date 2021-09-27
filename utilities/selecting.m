% [s] = selecting(cellbase, str)
%
% Takes a cellbase (meaning a cell array where the first row is all
% strings) and runs through each row; at each row, instantiates variables
% with the names of the first row giving them the value of the
% corresponding column in the current row; then runs eval(str), and puts
% the resulting value in the corresponding row of a vector s.
%
% Typically str will be some Boolean function.
%
% Example:
%    selecting(fof_list, 'strcmp(ratname, ''B068'') & lr>0.75')
%

function [my___s] = selecting(cellbase, my___str)

my___s = zeros(size(cellbase,1)-1,1);

for my___i=1:size(cellbase,1)-1,
	instantiate_row(cellbase(1,:), cellbase(my___i+1,:));
	my___s(my___i) = eval(my___str);
end;

return;


function instantiate_row(column_names, myrow)

for i=1:numel(column_names)
	assignin('caller', column_names{i}, myrow{i});
end;

