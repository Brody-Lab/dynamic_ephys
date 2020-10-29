function example_cell_psth(varargin)


p = inputParser;
addParameter(p,'cells',[ 16857 17784 18181 ]);
addParameter(p,'type','choice'); %evR, chrono
addParameter(p,'meta',0)
addParameter(p,'norm','none') %'none','onset','peak'
addParameter(p,'flip',0)
addParameter(p,'pause',0)
addParameter(p,'edges',[])


parse(p,varargin{:});
cell_to_plot = p.Results.cells;
type = p.Results.type;
meta = p.Results.meta;
norm = p.Results.norm;
flip = p.Results.flip;
do_pause = p.Results.pause;
edges = p.Results.edges;

if length(cell_to_plot) > 1 & ~meta
    for cc = 1:length(cell_to_plot)
        fprintf('working on cell %i (%i of %i)',...
            cell_to_plot(cc), cc, length(cell_to_plot))
        example_cell_psth('cells',cell_to_plot(cc),...
            'type',type,'meta',meta,'norm',norm,'flip',flip,...
            'edges',edges);
        if do_pause
            pause()
        end
    end
    return
end


dp = set_dyn_path;
%%
repack = 0;
ncells = length(cell_to_plot);

align_strs  = dyn_align_LUT;



min_t = 0;

gamma       = [];
T           = []; 
sides       = [];
poke_r      = [];
hit         = [];
evR         = [];
cin_fr      = [];
cout_fr     = [];
norm_f      = [];
end_state_s = [];



for cc = 1:ncells
    d       = dyn_cell_packager(cell_to_plot(cc),'repack',repack);
    
    if ~isfield(d.trials,'end_state_dur')
        d = dyn_cell_packager(cell_to_plot(cc),'repack',1);
    end
    
    coutind = find(ismember(d.align_strs,'cpokeout'));
    cinind  = find(ismember(d.align_strs,'stimstart-cout-mask'));
    stimind = strmatch('stimend',d.align_strs,'exact');
    if meta & p.Results.flip == 0
        flip(cc) = ~d.prefsideright{stimind};
    end
    T       = [T; d.trials.T]; %#ok<*AGROW>
    this_sides = 2.*d.trials.genEndState-1;
    this_poke_r = d.trials.rat_dir==1;
    this_gamma  = d.trials.gamma;
    if flip(cc)==1
        this_gamma = -this_gamma;
        this_sides = -this_sides;
        this_poke_r = ~this_poke_r;
    end
    gamma   = [gamma; this_gamma];

    sides   = [sides; this_sides];
    poke_r  = [poke_r; this_poke_r];
    hit     = [hit; d.trials.hit==1];
    evR     = [evR; d.trials.evidenceRatio(:)];
    end_state_s = [end_state_s; d.trials.end_state_dur(:)];

    switch norm
        case 'onset'
            norm_mult = d.norm_mean;
        case 'peak'
            error('doesn''t work yet')
        case 'none'
            norm_mult = 1;
    end
    norm_f  = [norm_f; norm_mult*ones(size(d.trials.T))];%;
    cin_fr  = [cin_fr; d.frate{cinind}];
    cout_fr = [cout_fr; d.frate{coutind}];
    
end
cin_t       = d.frate_t{cinind};
cout_t      = d.frate_t{coutind};


switch type
    case 'logR'
        logR = log(evR)./abs(gamma);
        if isempty(edges)
        edges = percentile(abs(logR),[30 70 ]);
        edges = sort([-Inf; -edges; 0; edges; Inf]);
        end
        [~, ~, bins] = histcounts(logR,edges);
        show_errors = 0;
        if meta
            error('not sure how to flip logR yet')
        end
    case 'choice'
        bins = poke_r + 1;
        show_errors = 1;
    case 'chrono'
        chrono = end_state_s .* sides;
        if isempty(edges)
                    nbins = 7;

        edges = linspace(-2,2,nbins);
        end
        [a, b, bins] = histcounts((chrono),edges);
        show_errors = 0;
        
end
nbins = length(unique(bins));

if meta
    cm = colormapRedBlue(ceil(nbins/2));
    cm(ceil(nbins/2)+1,:) = [];
else
cm = color_set(nbins);
end

err = hit == 0;

good = T > min_t;

fh = figure(1); clf

ax(1) = subplot(121);hold(ax(1),'on');
ax(2) = subplot(122);hold(ax(2),'on');

plot(ax(1),[ 0 0], [0 100],'k')
plot(ax(2),[ 0 0], [0 100],'k')

psths = [];

good_cint   = cin_t > -1 & cin_t < 1.5;
good_coutt  = cout_t > -1.5 & cout_t < .75;


for bb = 1:nbins
    this = good & bins == bb;
    this_color = cm(bb,:);

cin_hit_psth = nanmean(cin_fr(this  &  hit,good_cint)./norm_f(this&  hit));
cout_hit_psth = nanmean(cout_fr(this  &  hit,good_coutt)./norm_f(this&  hit));
psths = [psths ; cin_hit_psth cout_hit_psth];
if show_errors 
    err_color = this_color;
    %err_color = cm_(bb,:);
    cin_err_psth = nanmean(cin_fr(this  &  err,good_cint)./norm_f(this&  err));
    cout_err_psth = nanmean(cout_fr(this  & err,good_coutt)./norm_f(this&  err));
    psths = [psths ; cin_err_psth cout_err_psth];

    plot(ax(1),cin_t(good_cint),cin_err_psth,'--','color',err_color,'linewidth',1);
    plot(ax(2),cout_t(good_coutt),cout_err_psth,'--','color',err_color,'linewidth',1);

end

plot(ax(1),cin_t(good_cint),cin_hit_psth,'color',this_color,'linewidth',2);
plot(ax(2),cout_t(good_coutt),cout_hit_psth,'color',this_color,'linewidth',2);

end
linkaxes(ax,'y')

ylim(ax(1),([floor(min(psths(:))*10)/10 ceil(max(psths(:))*10)/10]))

xlim(ax(1),[cin_t(find(good_cint,1,'first')) cin_t(find(good_cint,1,'last'))])
xlim(ax(2),[cout_t(find(good_coutt,1,'first')) cout_t(find(good_coutt,1,'last'))])

%ax(2).YColor = 'w';
ax(1).TickDir = 'out';
ax(2).TickDir = 'out';

ylabel(ax(1), 'firing rate (spks/s)')
xlabel(ax(1), 'time from stim onset (s)') 
xlabel(ax(2), 'from movement (s)')
set(fh,'position',[5 5 6 3 ],'papersize',[5 3],'paperpositionmode','auto')

if ncells == 1
    group_name = num2str(cell_to_plot);
else
    group_name = 'group';
end
print(fh, fullfile(dp.psth_fig_dir, ['cell_' group_name ]),...
    '-dsvg', '-painters')

