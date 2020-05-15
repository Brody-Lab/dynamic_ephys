function [] = rat_variance(results, rats, just_good)

linear = [0.23 0.145 0.10 0.09];
step   = [0.35 0.195 0.116 0.1];
ratnames = {'H037','H066', 'H084','H129'};

good_cells = results.cell_dex;
celldex = zeros(length(good_cells),4);
celldex(:,1) =all(ismember(rats,'H037'),2);
celldex(:,2)= all(ismember(rats,'H066'),2);
celldex(:,3)= all(ismember(rats,'H084'),2);
celldex(:,4)= all(ismember(rats,'H129'),2);

%var = nanmean(results.var_explained,1);
var = nanmean(results.var_explained_static,1);
var(var < 0) = 0;
figure(1); clf; 
for dex=1:4
subplot(2,2,dex); hold on;
if just_good
    histogram(var(good_cells & celldex(:,dex)), 30,'facecolor','k','edgecolor','k');
    mid =nanmedian(var(good_cells & celldex(:,dex)));
    plot(mid, max([max(ylim) 15]), 'kv','markerfacecolor','k','markersize',10)
else
    histogram(var(logical(celldex(:,dex))), 30,'facecolor','k','edgecolor','k');
    mid = nanmedian(var(logical(celldex(:,dex))));
    plot(mid, max([max(ylim) 15]), 'kv','markerfacecolor','k','markersize',10)
    plot([mid mid], [0 max([max(ylim) 15])], 'k--')
    
end
xlim([0 0.5])

set(gca,'fontsize',12)
title([ratnames{dex} ' ' num2str(mid*100)],'fontsize',12)
ylabel('Count')
xlabel('variance explained')
pbaspect([1.25 1 1])
end

