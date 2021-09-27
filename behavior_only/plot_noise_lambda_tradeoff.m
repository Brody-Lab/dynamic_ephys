function [fh ax] = plot_noise_lambda_tradeoff(F,group);
dp = set_dyn_path;


fh = figure(1); clf;
ax = axes;
hold(ax ,'on');
load(fullfile(dp.model_fits_dir,'optimal_lambda'));
%yyaxis left
%ax = gca;
%ax.YColor = 'black';
msz = 12;
dp = set_dyn_path;
if ~strcmp(group, 'BING')
for i=1:length(C)
    plot(C{i}.noise, -C{i}.lambda,'.','color',[1 1 1].*.025,'markersize',.8*msz);
end
else
    plot([0 .5], [0 0], 'k-','linewidth',2);
end
%load('optimal_gaussian_lambda');
%for i=1:length(C)
%    plot(C{i}.noise, -C{i}.lambda,'marker','o','color',[.7 .7 .7],'markerfacecolor',[.7 .7 .7])
%end

for i=1:length(F)
    if isstruct(F{i})
        if F{i}.se(1) < 1
        plot(F{i}.n, -F{i}.final(1),'.-','color',[1 .55 .55],...
            'markerfacecolor',[1 .55 .55],'markersize',.8*msz)
        plot(F{i}.n+[F{i}.n_std, -F{i}.n_std], [-F{i}.final(1) -F{i}.final(1)],'-','color', dp.model_color)
        se = F{i}.se(1);
        plot([F{i}.n F{i}.n],-F{i}.final(1)+[se -se],'-','color', dp.model_color)
        end
    end
end
ylim([0 30])
xlim([0 .5])
if strcmp(group, 'BING')
ylim([-4 4])
end
ylabel('discounting rate (clks/s)');
xlabel('rat noise level (n)');

%yyaxis right
%ax = gca;
%ax.YColor = 'black';
%ylim([0 30])
%yticks([0 5 10 20 30])
%yticklabels({'\infty','0.2','0.1', '0.05', '0.03' })
%ylabel('$\tau ( = 1/\lambda$, sec)','Interpreter','latex','fontsize',20)

ax.TickDir = 'out';
axis('square');

fh.PaperUnits = 'inches';
pbaspect([1 1 1])
fh.Position = [3 3 3 3];
fh.PaperPosition = [0 0 3 3];
fh.PaperPositionMode = 'auto';
fh.PaperSize = [3 3];
print(fh,fullfile(dp.fig_dir, 'model_noise_v_lambda_tradeoff'),'-dsvg')



