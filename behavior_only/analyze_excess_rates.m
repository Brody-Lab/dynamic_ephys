function clicks = analyze_excess_rates(rat, S, p)

% log current step
disp(['Plotting Excess Click Rates for ' rat]);

% Analysis parameters
p.excess.win_dur    = 0.5;

if ~isfield(p.excess, 'OVERRIDE')
    p.excess.wStd       = 80;
    p.excess.ds         = 10;   
end 
%if p.haz == 0
%p.excess.wStd = 500;
%p.excess.ds = 5;
%end

p.excess.dt         = 0.0001;
p.excess.ww         = p.excess.wStd*4;
if p.haz == 0
    p.excess.mean       = 0;
else
    p.excess.mean       = 1;
end
p.excess.normalize  = 1;
p.excess.minT       = p.excess.win_dur + .05;
p.excess.window=exp(-.5.*(1:p.excess.ww).^2./p.excess.wStd^2)./sum(exp(-.5.*(1:p.excess.ww).^2./p.excess.wStd^2));
p.excess.full_steps = length(p.excess.dt:p.excess.dt:p.excess.minT);
p.excess.steps      = length(p.excess.dt:p.excess.dt*p.excess.ds:p.excess.minT);
p.compute_model     = 0;

% Compute behavior of linear model
p.optimal.compute   = 1;
p.optimal.lambda    = -4.1;
S = compute_linear_agent_dataset(rat,S,p);

% Call helper function that compute click rates
warning('off','all')
clicks = compute_excess_rates(S, p);


if 0
% make figure
if p.figvisible
    figure; hold on;
else
    figure('visible','off'); hold on;
end

if isstruct(clicks)

if p.excess.plotOptimal
    if length(clicks.RexcessOpt) ~=0
        shadedErrorBar(clicks.timepoint(1:50:end),clicks.RexcessOpt(1:50:end),clicks.RstdOpt(1:50:end), 'k',0)  
%    plot(clicks.timepoint(1:50:end), clicks.RexcessOpt(1:50:end),'k')
    end
    if length(clicks.LexcessOpt) ~=0
        shadedErrorBar(clicks.timepoint(1:50:end),clicks.LexcessOpt(1:50:end),clicks.LstdOpt(1:50:end), 'k',0)  
%    plot(clicks.timepoint(1:50:end), clicks.LexcessOpt(1:50:end),'k')
    end
end

% plot rat performance
if length(clicks.RexcessRat) ~=0
    shadedErrorBar(clicks.timepoint(1:50:end),clicks.RexcessRat(1:50:end),clicks.RstdRat(1:50:end), {'Color', p.nice_color{1}},0)  
end
if length(clicks.LexcessRat) ~=0
    shadedErrorBar(clicks.timepoint(1:50:end),clicks.LexcessRat(1:50:end),clicks.LstdRat(1:50:end), {'Color',p.nice_color{2}},0)  
end


if p.optimal.plotNoiseOptimal
    if length(clicks.RexcessNoise) ~=0
        shadedErrorBar(clicks.timepoint(1:50:end),clicks.RexcessNoise(1:50:end),clicks.RstdNoise(1:50:end), 'm',0)  
    end
    if length(clicks.LexcessNoise) ~=0
        shadedErrorBar(clicks.timepoint(1:50:end),clicks.LexcessNoise(1:50:end),clicks.LstdNoise(1:50:end), 'm',0)  
    end
end
if p.model.plot
    if length(clicks.RexcessModel) ~=0
        shadedErrorBar(clicks.timepoint(1:50:end),clicks.RexcessModel(1:50:end),clicks.RstdModel(1:50:end), p.model_color,0)  
    end
    if length(clicks.LexcessModel) ~=0
        shadedErrorBar(clicks.timepoint(1:50:end),clicks.LexcessModel(1:50:end),clicks.LstdModel(1:50:end), p.model_color,0)  
    end
end
if isfield(clicks,'timescaleRatFit') & p.excess.plotFit
    plot(clicks.timescaleRatFit)
end

axis square
legend(gca, 'off')
%title(rat, 'fontsize', p.singletitlesize)    
ylabel('Excess Rate', 'fontsize',p.mfontsize);
xlabel('Time from end of trial (s)', 'fontsize',p.mfontsize) 
xlim(clicks.xlim)
ylim([-p.excess.lim p.excess.lim])
if p.excess.normalize & p.excess.n_by_max
    ylim([-1.5 1.5])
end
set(gca, 'fontsize', p.mfontsize)
end

p.singleheight = 6;
p.singlewidth = 6;
% Save figure
if p.savefigure
     if strcmp(version('-release'), '2013b')
        fig = gcf;
        set(fig, 'PaperUnits', 'inches')
        set(fig, 'PaperPosition', [0 0 p.singlewidth p.singleheight])
        set(fig, 'PaperPositionMode', 'Manual')
    else
        fig = gcf;
        fig.PaperUnits = 'inches';
        fig.PaperPosition = [0 0 p.singlewidth p.singleheight];
        fig.PaperPositionMode = 'Manual';
        fig.PaperSize = [p.singlewidth p.singleheight];

    end 


%     if 7~=exist([p.figpath_root rat])
%         mkdir([p.figpath_root rat])
%     end
    fig_fn = fullfile(p.fig_dir,[ rat '_excess' p.hazStr]);
     
    print(fig_fn, '-dsvg')
    disp(['Figure saved to:' fig_fn]);

end

% Close figure
if ~p.figvisible
    close all;
end
end

% Save data for summary
haz = S.haz;

save_dir = fullfile(p.data_dir, rat);
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
save_fn = fullfile(save_dir, ['excess' p.hazStr]);
save(save_fn, 'clicks');
disp(['Data saved to:  ' save_fn]);

% Turn warnings back on
warning('on','all')
