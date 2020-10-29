p = set_dyn_path;
over = 1;
% load data
f = [];
d = [];
%%
for rr = 6:length(p.ratlist)
    ratname = p.ratlist{rr};
    fn = fullfile(p.behav_data_dir, [ratname '.mat']);
    fitfn = fullfile(p.model_fits_dir, ['fit_analytical_' ratname '.mat']);
    f{rr} = fit_rat_analytical(ratname,'results_dir', p.model_fits_dir,...
        'data_dir',p.behav_data_dir,'overwrite',over); %#ok<*SAGROW>
    %%
    d{rr} = load(fn,'data');
    data = d{rr}.data;
    %%
    assert(length(unique([data.Hazard]'))<5)
    %%
    endState = [data.genEndState]';
    switch_count = cellfun(@length, {data.genSwitchTimes})';
    no_switch = (switch_count==0);
    has_switch = switch_count >= 1;
    
    T = [data.T]';
    time_from_last = T;
    
    time_from_last(has_switch) = time_from_last(has_switch) -...
        cellfun(@(x) x(end), {data(has_switch).genSwitchTimes})';
    
    tbins = [0:.15: 1  2]
    
    cm = flipud(bone(length(tbins)));
    
    r = [data.pokedR]';
    h = [data.hit]';
    %% log odds psychometric
    
    gamma   = [data.gamma]';
    evR     = [data.evidenceRatio]';
    r       = [data.pokedR]';
    s_      = log(evR)./abs(gamma);
    ds      = .5;
    %sedges  = ceil(max(abs(s_))
    edges = [ -3:ds:3 ]
    sides = [data.correctAnswer]';
    medges = [ -3:.5:3 ]
    
    model_r = f{rr}.pr';
    
    model_h = (1 - model_r).*(~sides)+model_r.*sides;
    
    fh = figure(1); clf;
    set(fh,'position',[5 5 3 3]);
    
    ax = axes();
    plotPsychometric(s_, model_r, 'axHandle',ax,...
        'compute_fit', 0, 'plotfit',0, 'edges',medges,...
        'dataLineStyle','-','ploterrorbar',1,...
        'errorbar','gaussian',...
        'dataShaded',1, 'dataColor', p.model_color);
    
    plotPsychometric(s_, r, 'axHandle',ax,...
        'compute_fit', 0, 'plotfit',1,'edges',edges,...
        'dataLineStyle','.');
    
    ax.XGrid = 'off';
    ax.YGrid = 'off';
    ax.TickDir = 'out';
    ax.XTick = [-3 0 3];
    ax.XLim = [-3.1 3.1];
    ax.YLim = [-.02 1];

    ylabel('prob. go right');
    xlabel('log-odds supporting ''go right''');
    hl=legend('model','data','location','eastoutside');
    box(hl,'off');
    hl.Position = hl.Position + [.0 -.225 0 0];
    
    print(fh, fullfile(p.fig_dir, [ratname '_odds_psycho.svg']),'-dsvg','-painters')
    
    
    %% end state chrono
    %%
    dt = .2;
    fh = figure(3); clf
    ax = axes;
    
    
    ind = true(size(has_switch));
    tedges = [0:dt:max(time_from_last(ind))];
    set(fh,'position',[11 5 3 3],'papersize',[4 4])
    plotPsychometric(time_from_last(ind), model_h(ind),...
        'edges',tedges, 'axHandle',ax,'compute_fit',0,'plotfit',0,...
        'nbin',20,...
        'dataLineStyle','-','ploterrorbar',1,...
        'dataShaded',1, 'dataColor', p.model_color)
    
    plotPsychometric(time_from_last(ind), h(ind),...
        'edges',tedges, 'axHandle',ax,'compute_fit',0,'plotfit',1,...
        'dataLineStyle','.')
    
    xlabel('final state duration (s)')
    
    ax.XGrid = 'off';
    ax.YGrid = 'off';
    ax.TickDir = 'out'
    ylim([.48 1])
    xlim([-.1 .1] + tedges([1 end]));
    plot([-.1 -.1], [.48 .495], 'w', 'linewidth', 2.5)
    plot([-.1 -.005], [.48 .48], 'w', 'linewidth', 2.5)
    %
    % ax.XRuler.Axle.Visible = 'off'
    % ax.YRuler.Axle.Visible = 'off'
    % ax.YRuler.Axle.Visible = 'off'
    
    set(ax,'ytick',[.5:.1:1])
    ylabel('prob. correct')
    set(ax,'yticklabel',{'50' '' '' '' '' '100'})
    ylabel('prob. correct')
    set(ax,'yticklabel',{'.5' '' '' '' '' '1'})
    
    print(fh, fullfile(p.fig_dir, [ratname '_chrono.svg']),'-dsvg')
    
    
    %%
    dt = .2;
    fh = figure(3); clf
    ax = axes;
    ind = true(size(has_switch));
    tedges = [0:dt:max(time_from_last(ind))];
    set(fh,'position',[11 5 3 3],'papersize',[4 4])
    ctbins = [2 .8 0];
    cm = (bone(length(ctbins+1)))
    
    plotPsychometric(time_from_last(ind), model_h(ind),...
        'edges',tedges, 'axHandle',ax,'compute_fit',0,'plotfit',0,...
        'nbin',20,...
        'dataLineStyle','-','ploterrorbar',1,...
        'dataShaded',1, 'dataColor', p.model_color)
    
    for ii = 2:length(ctbins)
        ind = T <= ctbins(ii-1) & T > ctbins(ii);
        [~,res] = plotPsychometric(time_from_last(ind), h(ind),...
            'edges',tedges, 'axHandle',ax,'compute_fit',0,'plotfit',1,...
            'dataLineStyle','o','dataMarkerSize',6,'dataFaceColor',cm(ii-1,:))
    end
    xlabel('final state duration (s)')
    
    ax.XGrid = 'off';
    ax.YGrid = 'off';
    ax.TickDir = 'out'
    ylim([.48 1])
    xlim([-.1 .1] + tedges([1 end]));
    plot([-.1 -.1], [.48 .495], 'w', 'linewidth', 2.5)
    plot([-.1 -.005], [.48 .48], 'w', 'linewidth', 2.5)
    %
    % ax.XRuler.Axle.Visible = 'off'
    % ax.YRuler.Axle.Visible = 'off'
    % ax.YRuler.Axle.Visible = 'off'
    
    set(ax,'ytick',[.5:.1:1])
    ylabel('prob. correct')
    set(ax,'yticklabel',{'50' '' '' '' '' '100'})
    ylabel('prob. correct')
    set(ax,'yticklabel',{'.5' '' '' '' '' '1'})
    
    hl=legend('model','short trials','long trials','location','eastoutside')
    box(hl,'off')
    hl.Position = hl.Position + [.0 -.225 0 0]
    print(fh, fullfile(p.fig_dir, [ratname '_chrono_short_long.svg']),'-dsvg')
    
    %%
    dt = .2;

    ind = true(size(has_switch));
    tedges = [0:dt:max(time_from_last(ind))];
        fh = figure(3); clf
    ax = axes;
    set(fh,'position',[11 5 3 3],'papersize',[4 4])
    stbins = [0  1 2 Inf];
    cm = (bone(length(stbins)))
    
    plotPsychometric(time_from_last, model_h,...
        'edges',tedges, 'axHandle',ax,'compute_fit',0,'plotfit',0,...
        'nbin',20,...
        'dataLineStyle','-','ploterrorbar',1,...
        'dataShaded',1, 'dataColor', p.model_color)
    
    for ii = 2:length(stbins)
        
        ind = switch_count >= stbins(ii-1) & switch_count < stbins(ii);
        [~,res] = plotPsychometric(time_from_last(ind), h(ind),...
            'edges',tedges, 'axHandle',ax,'compute_fit',0,'plotfit',1,...
            'dataLineStyle','o','dataMarkerSize',6,'dataFaceColor',cm(ii-1,:))
        
    end
    xlabel('final state duration (s)')
    
    ax.XGrid = 'off';
    ax.YGrid = 'off';
    ax.TickDir = 'out'
    ylim([.48 1])
    xlim([-.1 .1] + tedges([1 end]));
    plot([-.1 -.1], [.48 .495], 'w', 'linewidth', 2.5)
    plot([-.1 -.005], [.48 .48], 'w', 'linewidth', 2.5)
    %
    % ax.XRuler.Axle.Visible = 'off'
    % ax.YRuler.Axle.Visible = 'off'
    % ax.YRuler.Axle.Visible = 'off'
    
    set(ax,'ytick',[.5:.1:1])
    ylabel('prob. correct')
    set(ax,'yticklabel',{'50' '' '' '' '' '100'})
    ylabel('prob. correct')
    set(ax,'yticklabel',{'.5' '' '' '' '' '1'})
    
    hl=legend('model','0 switches','1','many','location','eastoutside')
    box(hl,'off')
    hl.Position = hl.Position + [.0 -.225 0 0]
    print(fh, fullfile(p.fig_dir, [ratname '_chrono_nswitches.svg']),'-dsvg',...
        '-painters')
    
    
    
    
    %%
    fh = figure(5);clf
    set(fh,'position',[11 1 3 3],'papersize',[4 4])
    
    ax = axes
    s_ = [data.Delta]';
    side_r = sides == 1;
    xmax = round(max(abs(s_))/10)*10
    dedges = linspace(-xmax,xmax,9)
    plotPsychometric(s_, side_r, 'axHandle',ax,...
        'compute_fit', 0, 'plotfit',0,...
        'dataColor', 'r','edges',dedges)
    plotPsychometric(s_, r, 'axHandle',ax,...
        'compute_fit', 0, 'plotfit',0,...
        'dataColor', 'k','edges',dedges)
    % plotPsychometric(s_(has_switch), side_r(has_switch), 'axHandle',ax,...
    %     'compute_fit', 0, 'plotfit',0,...
    %     'dataColor', one_switch_color,'edges',dedges)
    xlabel('click difference (r-l)')
    ylabel('prob. go right')
    hl=legend('true side','rat choice','no','location','east')
    title(hl,'')
    box(hl,'off')
    ax.YGrid = 'off';
    ax.XGrid = 'off';
    hl.Position = hl.Position + [.05 -.225 0 0]
    %%
end


%%
% population chrono
fh = figure(4); clf
ax = axes;
set(fh,'position',[11 5 3 3],'papersize',[4 4])
pop_tfl = []; pop_h = [];
for dd = 1:length(d)
    data = d{dd}.data;
        switch_count = cellfun(@length, {data.genSwitchTimes})';
    no_switch = (switch_count==0);
    has_switch = switch_count >= 1;
    T = [data.T]';
    time_from_last = T;
    r = [data.pokedR]';
    h = [data.hit]';
    
    time_from_last(has_switch) = time_from_last(has_switch) -...
        cellf   un(@(x) x(end), {data(has_switch).genSwitchTimes})';
    plotPsychometric(time_from_last, h,...
        'edges',tedges, 'axHandle',ax,'compute_fit',0,'plotfit',1,...
        'ploterrorbar',0,'dataLineStyle','-','dataLineWidth',1,...
        'dataColor',[1 1 1].*.5)
    
    pop_tfl{dd} = time_from_last;
    pop_h{dd}   = h;
end
pop_tfl = vertcat(pop_tfl{:});
pop_h = vertcat(pop_h{:});
plotPsychometric(pop_tfl, pop_h,...
        'edges',tedges, 'axHandle',ax,'compute_fit',0,'plotfit',1,...
        'ploterrorbar',0,'dataLineStyle','-','dataLineWidth',2,...
        'dataColor',[1 1 1].*0)
xlabel('final state duration (s)')
    
    ax.XGrid = 'off';
    ax.YGrid = 'off';
    ax.TickDir = 'out'
    ylim([.48 1])
    xlim([-.1 .1] + tedges([1 end]));
    plot([-.1 -.1], [.48 .495], 'w', 'linewidth', 2.5)
    plot([-.1 -.005], [.48 .48], 'w', 'linewidth', 2.5)
    %
    % ax.XRuler.Axle.Visible = 'off'
    % ax.YRuler.Axle.Visible = 'off'
    % ax.YRuler.Axle.Visible = 'off'
    
    set(ax,'ytick',[.5:.1:1])
    ylabel('prob. correct')
    set(ax,'yticklabel',{'50' '' '' '' '' '100'})
    ylabel('prob. correct')
    set(ax,'yticklabel',{'.5' '' '' '' '' '1'})
    
    %hl=legend('model','0 switches','1','many','location','eastoutside')
    %box(hl,'off')
    %hl.Position = hl.Position + [.0 -.225 0 0]
    print(fh, fullfile(p.fig_dir, ['population_chrono.svg']),'-dsvg',...
        '-painters')
    
    