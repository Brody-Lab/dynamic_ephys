function [results] = decompose_tuning_curves(results,cellids)

if nargin == 1
    cellids = results.cellid;
end
if ~isfield(results, 'tuning_cell')
    results.tuning_cell = NaN(size(results.fga_cell_residual));%svd approx
    results.tuning_fr = NaN(size(results.frm_time));%vector of fr modulation
    results.tuning_fa = NaN(size(results.fga_ta_cell));%tuning curve 
    results.rank1_variance = NaN(size(results.cellid));
end
for i=1:length(cellids)
    results = decompose_tuning_curves_single(results, cellids(i));
end
end

function [results] = decompose_tuning_curves_single(results, cellid)
    celldex = find(results.cellid == cellid);
    if isempty(celldex)
        keyboard
        error('Couldnt find cell in cellid list')
    end
    fga = results.fga_cell_residual(:,:,celldex);
    [u,s,v] = svd(fga);

    % compute Rank 1 approx
    s2 = s;
    s2(2:end) = 0;
    results.tuning_cell(:,:,celldex) = u*s2*v';
    results.tuning_fr(:,celldex) = u(:,1);
    results.tuning_fa(:,celldex) = v(:,1);
    s_squared = diag(s).^2;
    results.rank1_variance(celldex) = s_squared(1)./sum(s_squared);
%    results.rank1_variance(celldex) = s(1)./sum(diag(s)); 
end
