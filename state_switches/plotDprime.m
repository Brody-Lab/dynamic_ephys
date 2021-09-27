function plotDprime(cellid, which_switch, savefig)

if nargin < 2
    which_switch = 'model';
end


%%% look at a neuron of interest
figure
res =  compute_switch_triggered_average(cellid,'post',2,'which_switch',which_switch, 'n_shuffles', 1000,'save_file',1,'mask_other_switch',1);
plot_sta_dprime(res, gcf)

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 6];
fig.PaperPositionMode = 'Manual';
fig.PaperSize = [6 6];
if nargin > 2 && savefig
    print(gcf,  ['../figures/STA/' num2str(cellid) '_' which_switch '_STA_dprime'],'-dsvg')
end



