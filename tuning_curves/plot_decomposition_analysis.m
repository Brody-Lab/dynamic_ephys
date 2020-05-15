function [] = plot_decomposition_analysis(results)
    figure;
    subplot(2,2,1); hold on;
    plot(1:length(results.cellid), flipud(sort(results.rank1_variance).*100), 'ko')
    ylabel('Rank 1 variance')
    xlabel('Sorted #')
    set(gca, 'fontsize',16)
    ylim([0, 100])
    xlim([0, length(results.cellid)+1])
%    subplot(2,2,2); hold on;
%    plot(1:length(results.cellid), fliplr(sort(results.rank1_variance_squared).*100), 'ko')
%    ylabel('Rank 1 variance')
%    xlabel('Sorted #')
%    set(gca, 'fontsize',16)
%    ylim([0, 100])
%    xlim([0, length(results.cellid)+1])

    subplot(2,2,2); hold on;
    plot(results.fr_modulation, results.rank1_variance.*100, 'ko')
    ylabel('Rank 1 variance')
    xlabel('FR. Modulation (Hz)')
    set(gca, 'fontsize',16)
    ylim([0, 100])
    subplot(2,2,3); hold on;
    plot(results.svd_slope_cell, results.rank1_variance.*100, 'ko')
    ylabel('Rank 1 variance')
    xlabel('Tuning curve slope')
    set(gca, 'fontsize',16)
    ylim([0, 100])
end
