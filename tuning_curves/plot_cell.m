function [] = plot_cell(results, varargin)

celldex = 1;
savefig = 0;
plot_sigmoid = 1;
plot_nice = 1;
res = 5;
if length(varargin) > 0 
    celldex = varargin{1};
    if numel(celldex) > 1
        % more than one cell, recursively plot them
        plot_cell(results, celldex(2:end))
        celldex = celldex(1);
    end
    if (celldex > length(results.cell_dex) ) | (celldex < 0)
        % was given cellid, rather than cell index
        newdex = find(results.cellid == celldex);
        if isempty(newdex)
            error('Bad cell index, neither logical nor cell id')
        else
            celldex = newdex;
        end
    end
    if nargin > 2
        savefig = varargin{2};
        if savefig        
            disp('Saving figures')
        end
    end
end
cellid = results.cellid(celldex);
if nargin > 3
    res = varargin{3};
end


if ndims(results.fga_cell) == 2
    fga_cell            = results.fga_cell;
    fga_cell_residual   = results.fga_cell_residual;
    fga_cell_nresidual  = results.fga_cell_nresidual;
    fga_ta_cell         = results.fga_ta_cell;
    fga_nta_cell        = results.fga_nta_cell;
    fga_std_ta_cell     = results.fga_std_ta_cell;
    fga_std_nta_cell    = results.fga_std_nta_cell;
    fga_ta_unorm_cell   = results.fga_ta_unorm_cell;
    fga_aa_cell         = results.fga_aa_cell;
    betas               = results.betas_cell;
    sigma               = results.sigmas_cell;
    nbetas              = results.nbetas_cell;
    nsigma              = results.nsigmas_cell;
    flipcond            = results.flipdex;
    frm_time            = results.frm_time;
    %var_expl            = results.var_explained;
    tuning_cell         = results.tuning_cell;
    tuning_fr           = results.tuning_fr;
    tuning_fa           = results.tuning_fa;
    svd_betas           = results.svd_betas_cell;
    svd_sigma           = results.svd_sigmas_cell;
else
    fga_cell            = results.fga_cell(:,:,celldex);
    fga_cell_residual   = results.fga_cell_residual(:,:,celldex);
    fga_cell_nresidual  = results.fga_cell_nresidual(:,:,celldex);
    fga_ta_cell         = results.fga_ta_cell(:,celldex);
    fga_nta_cell        = results.fga_nta_cell(:,celldex);
    fga_std_ta_cell     = results.fga_std_ta_cell(:,celldex);
    fga_std_nta_cell    = results.fga_std_nta_cell(:,celldex);
    fga_ta_unorm_cell   = results.fga_ta_unorm_cell(:,celldex);
    fga_aa_cell         = results.fga_aa_cell(:,celldex);
    betas               = results.betas_cell(celldex,:);
    sigma               = results.sigmas_cell(celldex,:);
    nbetas              = results.nbetas_cell(celldex,:);
    nsigma              = results.nsigmas_cell(celldex,:);
    flipcond            = results.flipdex(celldex);
    frm_time            = results.frm_time(:,celldex);
    %var_expl            = results.var_explained(:,celldex);
    tuning_cell         = results.tuning_cell(:,:,celldex);
    tuning_fr           = results.tuning_fr(:,celldex);
    tuning_fa           = results.tuning_fa(:,celldex);
    svd_betas           = results.svd_betas_cell(celldex,:);
    svd_sigma           = results.svd_sigmas_cell(celldex,:);
end

if ~flipcond
    fga_cell            = flipdim(fga_cell,2);
    fga_cell_residual   = flipdim(fga_cell_residual,2);
    fga_cell_nresidual  = flipdim(fga_cell_nresidual,2);
    fga_ta_cell         = flipdim(fga_ta_cell,1);
    fga_nta_cell        = flipdim(fga_nta_cell,1);
    fga_ta_unorm_cell   = flipdim(fga_ta_unorm_cell,1); 
    tuning_cell         = flipdim(tuning_cell,2);
    tuning_fa           = flipdim(tuning_fa,1);

end
% plot fga_ta_cell, with sigmoid for both cases
n_dv_bins = numel(fga_ta_cell);
%clrs            = rainbow_colors(n_dv_bins);
clrs = color_set(n_dv_bins);
markersize = 2;

figure(1); clf;
%%%% Plot r(a,t) 
subplot(3,4,1);hold on; set(gca, 'fontsize',12);
for i=1:res:n_dv_bins
    if res > 1
        plot(results.t0s, movmean(fga_cell(:,i),res), '-', 'color', clrs(i,:), 'Markerfacecolor', clrs(i,:), 'MarkerSize', markersize);

    else
        plot(results.t0s, fga_cell(:,i), '-', 'color', clrs(i,:), 'Markerfacecolor', clrs(i,:), 'MarkerSize', markersize);
    end
end
ylabel('Firing Rate (Hz)')
xlabel('Time (s)')
title(['Cell ' num2str(cellid) ' r(a,t)'])
%pbaspect([1.25 1 1])
ylim([0 inf])
ylims1 = ylim;

%%%% Plot Delta r(a,t) 
subplot(3,4,2);hold on; set(gca, 'fontsize',12);
for i=1:res:n_dv_bins
    if res > 1
    plot(results.t0s, movmean(fga_cell_residual(:,i),res), '-', 'color', clrs(i,:), 'Markerfacecolor', clrs(i,:), 'MarkerSize', markersize);
    else
    plot(results.t0s, fga_cell_residual(:,i), '-', 'color', clrs(i,:), 'Markerfacecolor', clrs(i,:), 'MarkerSize', markersize);
    end
end
ylabel('\Delta FR (Hz)')
xlabel('Time (s)')
title('\Delta r(a,t)')
%pbaspect([1.25 1 1])
ylims = ylim;
ylim([-max(abs(ylims)), +max(abs(ylims))])

%%%% Plot Normalized Delta r(a,t) 
subplot(3,4,4);hold on; set(gca, 'fontsize',12);
if plot_nice
    for i=1:n_dv_bins
        fga_cell_nresidual(:,i) = movmean(fga_cell_nresidual(:,i),res);
    end
    for i=1:length(results.t0s)
        fga_cell_nresidual(i,:) = fga_cell_nresidual(i,:) - min(fga_cell_nresidual(i,:));
        fga_cell_nresidual(i,:) = fga_cell_nresidual(i,:)./max(fga_cell_nresidual(i,:));
    end
end
for i=1:res:n_dv_bins
%    if res > 1
%    plot(results.t0s, movmean(fga_cell_nresidual(:,i),res+2), '-', 'color', clrs(i,:), 'Markerfacecolor', clrs(i,:), 'MarkerSize', markersize);
%    else
    plot(results.t0s(1:res:end), fga_cell_nresidual(1:res:end,i), '-', 'color', clrs(i,:), 'Markerfacecolor', clrs(i,:), 'MarkerSize', markersize);
%    end
end
ylabel('Norm. FR')
xlabel('Time (s)')
title('\Delta r(a,t)')
%pbaspect([1.25 1 1])
%ylims = ylim;
%ylim([-max(abs(ylims)), +max(abs(ylims))])


%plot(results.t0s, fga_aa_cell,'k','linewidth',2)
%ylabel('Firing Rate (Hz)')
%xlabel('Time (s)')
%title('Mean fr | t ')
%pbaspect([1.25 1 1])
%ylim(ylims1)

%%%% Plot Rank 1 Delta r(a,t) 
subplot(3,4,3);hold on; set(gca, 'fontsize',12);
for i=1:res:n_dv_bins
    if res > 1
    plot(results.t0s, movmean(tuning_cell(:,i),res), '-', 'color', clrs(i,:), 'Markerfacecolor', clrs(i,:), 'MarkerSize', markersize);
    else
    plot(results.t0s, tuning_cell(:,i), '-', 'color', clrs(i,:), 'Markerfacecolor', clrs(i,:), 'MarkerSize', markersize);
    end
end
ylabel('\Delta FR (Hz)')
xlabel('Time (s)')
title('\Delta r(a,t)-r1')
%pbaspect([1.25 1 1])
ylims = ylim;
ylim([-max(abs(ylims)), +max(abs(ylims))])


%%%% Plot tuning curve from avg method?
subplot(3,4,5);hold on; set(gca, 'fontsize',12);
plot(results.dv_axis,fga_ta_cell,'k','linewidth',2)
h = gca;
for i=1:res:numel(results.dv_axis)
    hold on;
    eh = errorplot2(h,results.dv_axis(i), fga_ta_cell(i), fga_std_ta_cell(i), 'Marker', '', 'Color', clrs(i,:));
end;
xfine = results.dv_axis(1):0.1:results.dv_axis(end);
if ~isnan(betas)
ypred = dyn_sig(betas,xfine);
if ~flipcond
    ypred = flipdim(ypred,2);
end
if plot_sigmoid
    plot(xfine, ypred, 'm','linewidth',2)
end
end
title('Avg. r(a)')
xlabel('Accumulated Value')
ylabel('Norm. FR')
%pbaspect([1.25 1 1])
ylim([-0.1 1.1])

%%%% Plot un-normalized tuning curve 
subplot(3,4,6);hold on; set(gca, 'fontsize',12);
plot(results.dv_axis,fga_ta_unorm_cell,'k','linewidth',2)
%fr_mod = 0.75;
%plot([results.dv_axis(1) results.dv_axis(end)], [fr_mod/2 fr_mod/2], 'r--');
%plot([results.dv_axis(1) results.dv_axis(end)], [-fr_mod/2 -fr_mod/2], 'r--');
plot([results.dv_axis(1) results.dv_axis(end)], [0 0], 'k--');
title('raw r(a)')
xlabel('Accumulated Value')
ylabel('\Delta FR (Hz)')
%pbaspect([1.25 1 1])
ylims = ylim;
ylim([-max(abs(ylims)), +max(abs(ylims))])
%
%%%% Plot normalized tuning curve 
subplot(3,4,7); hold on; set(gca, 'fontsize',12)
%title(['Cell ' num2str(cellid)])
plot(results.dv_axis,fga_nta_cell-0.5,'k','linewidth',2)
h = gca;
for i=1:res:numel(results.dv_axis)
    hold on;
    eh = errorplot2(h,results.dv_axis(i), fga_nta_cell(i)-0.5, fga_std_nta_cell(i), 'Marker', '', 'Color', clrs(i,:));
end;
xfine = results.dv_axis(1):0.1:results.dv_axis(end);
if ~isnan(nbetas)
ypred = dyn_sig(nbetas,xfine);
if ~flipcond
    ypred = flipdim(ypred,2);
end
if plot_sigmoid
    plot(xfine, ypred-0.5, 'm','linewidth',2)
end
end
title('Norm. r(a)')
xlabel('Accumulation Value (a)')
ylabel('Norm. FR')
%pbaspect([1.25 1 1])
%ylim([-0.1 1.1])
ylim([-0.6 0.6])
set(h, 'Ytick', [-0.5 0 0.5], 'yticklabel',{'-0.5','0','0.5'})

%%%% Plot SVD tuning curve
subplot(3,4,8); hold on; set(gca, 'fontsize',12)
plot([results.dv_axis(1) results.dv_axis(end)], [0 0], 'k--');
%title(['Cell ' num2str(cellid)])
plot(results.dv_axis,tuning_fa,'k','linewidth',2)
h = gca;

xfine = results.dv_axis(1):0.1:results.dv_axis(end);
if ~isnan(svd_betas)
ypred = dyn_sig(svd_betas,xfine);
if ~flipcond
    ypred = flipdim(ypred,2);
end
if plot_sigmoid
    plot(xfine, ypred+min(tuning_fa), 'm','linewidth',2)
end
end
title('SVD r(a)')
xlabel('Accumulation Value (a)')
ylabel('FR')
%pbaspect([1.25 1 1])
%ylim([-0.1 1.1])
ylim([-0.6 0.6])
set(h, 'Ytick', [-0.5 0 0.5], 'yticklabel',{'-0.5','0','0.5'})

%%%% Plot residual fr
subplot(3,4,9); hold on;
plot(results.t0s, fga_aa_cell, 'k','linewidth',2)
plot(results.t0s, zeros(size(results.t0s)),'Color',[.75 .75 .75],'linewidth',1);
title('Avg. FR')
ylabel('Avg. FR.(Hz)')
xlabel('Time (s)');
set(gca, 'fontsize',12);
%pbaspect([1.25 1 1])

%%%% Plot firing rate modulation
subplot(3,4,10); hold on;
plot(results.t0s, frm_time, 'k','linewidth',2)
plot(results.t0s, zeros(size(results.t0s)),'Color',[.75 .75 .75],'linewidth',1);
title('FR. Mod')
ylabel('FR. Mod. (Hz)')
xlabel('Time (s)');
set(gca, 'fontsize',12);
%pbaspect([1.25 1 1])

%%%% Plot firing rate modulation
subplot(3,4,11); hold on;
plot(results.t0s, tuning_fr, 'k','linewidth',2)
plot(results.t0s, zeros(size(results.t0s)),'Color',[.75 .75 .75],'linewidth',1);
title('SVD FR. Mod')
ylabel('FR mod. (Hz)')
xlabel('Time (s)');
set(gca, 'fontsize',12);
%pbaspect([1.25 1 1])

subplot(3,4,12); hold on;
%plot(results.t0s, var_expl.*100, 'k','linewidth',2)
plot(results.t0s, zeros(size(results.t0s)),'Color',[.75 .75 .75],'linewidth',1);
ylabel('VE %')
xlabel('Time (s)');
set(gca, 'fontsize',12);
%pbaspect([1.25 1 1])




if savefig
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 10 6];
    fig.PaperPositionMode = 'Manual';
    fig.PaperSize = [10 6];
    print(gcf,  ['../../figures/tuning/' num2str(cellid) ],'-dsvg')
end

if 0
[u,s,v] = svd(fga_cell_residual);
figure; plot(results.t0s, movmean(u(:,1).*s(1,1).*0.12,10),'k','linewidth',2)
ylabel('f.r modulation (Hz)')
xlabel('Time (s)')
set(gca, 'fontsize',12)
title('m(t)')
pbaspect([1.25 1 1])
ylim([-5 5])
if savefig
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 4 4];
    fig.PaperPositionMode = 'Manual';
    fig.PaperSize = [4 4];
    print(gcf,  ['../../figures/tuning/m_' num2str(cellid) ],'-dsvg')
end
end






