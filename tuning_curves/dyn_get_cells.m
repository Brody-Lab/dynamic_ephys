function [singles all_cells] = dyn_get_cells()

ratnames = {'H037', 'H066', 'H084', 'H129'};
all_cells =[];
single = [];
for rr = 1:length(ratnames)
    [single{rr}, all_cells{rr}] = bdata('select single, cellid from cells where ratname="{S}"',ratnames{rr});
end
single = vertcat(single{:});
all_cells = vertcat(all_cells{:});
singles = all_cells(single==1);