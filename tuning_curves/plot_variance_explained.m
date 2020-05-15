function [] = plot_variance_explained(results)

figure;
vals = results.var_explained_static(results.cell_dex);
vals(vals < 0) = [];
histogram(vals.*100, 30, 'facecolor','k','edgecolor','k')
ylabel('count')
xlabel('Variance Explained %')
set(gca, 'fontsize',12)

figure;
vals = results.var_explained(results.cell_dex);
vals(vals < 0) = [];
histogram(vals.*100, 30, 'facecolor','k','edgecolor','k')
ylabel('count')
xlabel('Variance Explained %')
set(gca, 'fontsize',12)


