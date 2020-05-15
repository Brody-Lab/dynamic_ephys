function [] = rat_histogram(results, rats, just_good)

linear = [0.23 0.145 0.10 0.09];
step   = [0.35 0.195 0.116 0.1];
ratnames = {'H037','H066', 'H084','H129'};

good_cells = results.cell_dex;
celldex = zeros(length(good_cells),4);
celldex(:,1) =all(ismember(rats,'H037'),2);
celldex(:,2)= all(ismember(rats,'H066'),2);
celldex(:,3)= all(ismember(rats,'H084'),2);
celldex(:,4)= all(ismember(rats,'H129'),2);

figure(1); clf; 
for dex=1:4
subplot(2,2,dex); hold on;
if just_good
    histogram(results.slope_cell(good_cells & celldex(:,dex)), 30,'facecolor','k','edgecolor','k');
    plot(nanmedian(results.slope_cell(good_cells & celldex(:,dex))), max([max(ylim) 15]), 'kv','markerfacecolor','k','markersize',10)
else
    histogram(results.slope_cell(logical(celldex(:,dex))), 30,'facecolor','k','edgecolor','k');
    mid = nanmedian(results.slope_cell(logical(celldex(:,dex))));
    plot(mid, max([max(ylim) 15]), 'kv','markerfacecolor','k','markersize',10)
    plot([mid mid], [0 max([max(ylim) 15])], 'k--')
    
end
plot(linear(dex), max([max(ylim) 15]), 'rv','markerfacecolor','r','markersize',10)
plot(step(dex), max([max(ylim) 15]), 'bv','markerfacecolor','b','markersize',10)
xlim([0 inf])

set(gca,'fontsize',12)
title(ratnames{dex},'fontsize',12)
ylabel('Count')
xlabel('svd slope')
pbaspect([1.25 1 1])
end

