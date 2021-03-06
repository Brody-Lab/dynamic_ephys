function [fh ax] = example_cell_psth(varargin)
% function [fh ax] = example_cell_psth(varargin) 
% example usage 
% [fh, ax] = example_cell_psth('cells', cellid,...
%     'cintrange', xlim_on, 'couttrange', xlim_off, ...
%     'coutstr', 'stimend-nomask', 'fig_num', 2)

p = inputParser;
addParameter(p,'cells',[ 16857 17784 18181 ]);
addParameter(p,'type','choice'); %evR, chrono
addParameter(p,'separate_hits',1);
addParameter(p,'meta',0)
addParameter(p,'norm','none') %'none','onset','peak'
addParameter(p,'flip',0)
addParameter(p,'pause',0)
addParameter(p,'top_color','')
addParameter(p,'bot_color','')
addParameter(p,'edges',[])
addParameter(p,'min_t',0)
addParameter(p,'fig_num',1)
addParameter(p,'repack',0)
addParameter(p,'cintrange',[-1 1.5])
addParameter(p,'couttrange',[-1.5 .75])
addParameter(p,'fpos',[5 5 6 3])
addParameter(p, 'cinstr', 'stimstart-cout-mask');
addParameter(p, 'coutstr', 'cpokeout');
addParameter(p, 'ploterrorbar', 0);
addParameter(p, 'errorbarfun', @nansem);
addParameter(p, 'dyn_path',  []);

parse(p,varargin{:});
cintrange = p.Results.cintrange;
couttrange = p.Results.couttrange;
cinstr = p.Results.cinstr;
coutstr = p.Results.coutstr;

cell_to_plot = p.Results.cells;
type = p.Results.type;
meta = p.Results.meta;
norm = p.Results.norm;
flip = p.Results.flip;
fpos = p.Results.fpos;
do_pause = p.Results.pause;
edges = p.Results.edges;
separate_hits = p.Results.separate_hits;
fig_num     = p.Results.fig_num;
ploterrorbar = p.Results.ploterrorbar;
errorbarfun = p.Results.errorbarfun;
dp          = p.Results.dyn_path;
repack      = p.Results.repack;


if isempty(dp)
    dp = set_dyn_path;
end

if length(cell_to_plot) > 1 & ~meta
    for cc = 1:length(cell_to_plot)
        fprintf('working on cell %i (%i of %i)',...
            cell_to_plot(cc), cc, length(cell_to_plot))
        example_cell_psth('cells',cell_to_plot(cc),...
            'type',type,'meta',meta,'norm',norm,'flip',flip,...
            'edges',edges,'ploterrorbar',ploterrorbar,...
            'dyn_path', dp, 'repack', repack);
        if do_pause
            pause();
            disp(cell_to_plot)
        end
    end
    return
end



%%

ncells = length(cell_to_plot);

align_strs  = dyn_align_LUT;

min_t = p.Results.min_t;

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

fn          = sprintf('meta_psth_%s.mat',type);
savename    = fullfile(dp.spikes_dir, fn);
if meta & exist(savename, 'file')
    load(savename)
else
    
    for cc = 1:ncells
        fprintf('cell %i (%i/%i)...\n',     cell_to_plot(cc), cc, ncells);
        d       = dyn_cell_packager(cell_to_plot(cc),'dyn_path',dp,'repack',repack);
        
        if ~isfield(d.trials,'end_state_dur')
            d = dyn_cell_packager(cell_to_plot(cc),'repack',1);
        end
        
        end_tind = find(ismember(d.align_strs,coutstr));
        start_tind  = find(ismember(d.align_strs,cinstr));
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
        cin_fr  = [cin_fr; d.frate{start_tind}];
        cout_fr = [cout_fr; d.frate{end_tind}];
        
    end
    cin_t       = d.frate_t{start_tind};
    cout_t      = d.frate_t{end_tind};
    save(savename, 'gamma', 'T', 'sides', 'poke_r', 'hit', 'evR', ...
        'cin_t', 'cin_fr', 'cout_t', 'cout_fr', 'norm_f', 'end_state_s');
    
end

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
        cblab = 'evR';
    case 'choice'
        bins = poke_r + 1;
        show_errors = 1;
        cblab = 'choice';
    case 'chrono'
        chrono = end_state_s .* sides;
        if isempty(edges)
            nbins = 7;
            
            edges = linspace(-2,2,nbins);
        end
        [a, b, bins] = histcounts((chrono),edges);
        show_errors = 0;
        cblab = {'final state' 'duration (s)'};
        
end
nbins = length(unique(bins));

if meta
    %%
    top = [.65 .2 .75];
    top = [.75 .4 .75];
    if ~isempty(p.Results.top_color)
        top = p.Results.top_color;
    end
    bot = [1 1 1] - top;
    if ~isempty(p.Results.bot_color)
        bot = p.Results.bot_color;
    end
    bot = [.75 .75 .25];
   
    cm_bot = colormapLinear(bot, ceil(nbins/2));
    cm_top = colormapLinear(top, ceil(nbins/2));
    cm = [cm_bot(end:-1:2,:); cm_top(2:end,:)];
    
else
    cm = color_set(nbins);
end

err = hit == 0;

good = T > min_t;

fh = figure(fig_num); clf
set(fh,'position',fpos,'papersize',fpos([3 4]),'paperpositionmode','auto')

ax(1) = subplot(121);hold(ax(1),'on');
ax(2) = subplot(122);hold(ax(2),'on');

plot(ax(1),[ 0 0], [0 100],'k')
plot(ax(2),[ 0 0], [0 100],'k')

psths       = [];
 
good_cint   = cin_t > cintrange(1) & cin_t < cintrange(2);
good_coutt  = cout_t > couttrange(1) & cout_t < couttrange(2);


for bb = 1:nbins
    this = good & bins == bb;
    this_color = cm(bb,:);
    if separate_hits
        trials = this & hit;
    else
        trials = this;
    end
    
    cin_hit_psth = nanmean(cin_fr(trials,good_cint)./norm_f(trials));
    cin_hit_psth_sem = errorbarfun(cin_fr(trials,good_cint)./norm_f(trials));
    cout_hit_psth = nanmean(cout_fr(trials,good_coutt)./norm_f(trials));
    cout_hit_psth_sem = errorbarfun(cout_fr(trials,good_coutt)./norm_f(trials));
    psths = [psths ; cin_hit_psth cout_hit_psth];
    if show_errors
        err_color = this_color;
        %err_color = cm_(bb,:);
        cin_err_psth = nanmean(cin_fr(this  &  err,good_cint)./norm_f(this&  err));
        cin_err_psth_sem = errorbarfun(cin_fr(this  &  err,good_cint)./norm_f(this&  err));
        cout_err_psth = nanmean(cout_fr(this  & err,good_coutt)./norm_f(this&  err));
        cout_err_psth_sem = errorbarfun(cout_fr(this  & err,good_coutt)./norm_f(this&  err));
        psths = [psths ; cin_err_psth cout_err_psth];
        
%         
%         plot(ax(1),cin_t(good_cint),cin_err_psth,'--','color',err_color,'linewidth',1.5);
%         plot(ax(2),cout_t(good_coutt),cout_err_psth,'--','color',err_color,'linewidth',1.5);
%         
            if ploterrorbar
                %%
                hsvcolor = rgb2hsv(this_color);
                this_color_light = hsv2rgb([hsvcolor(1) .35 .9]);
                set(fh,'currentaxes',ax(1))
                shadedErrorBar(cin_t(good_cint),cin_err_psth,cin_err_psth_sem,...
                    {'linestyle','--','color',this_color_light,'linewidth',1.5,'parent', ax(1)},1);
                set(fh,'currentaxes',ax(2))
                shadedErrorBar(cout_t(good_coutt),cout_err_psth, cout_err_psth_sem,...
                    {'linestyle','--','color',this_color_light,'linewidth',1.5,'parent', ax(2)},1);
            else
                plot(ax(1),cin_t(good_cint),cin_err_psth,'--','color',err_color,'linewidth',1);
                plot(ax(2),cout_t(good_coutt),cout_err_psth,'--','color',err_color,'linewidth',1);
            end
        
    end
    
    if ploterrorbar
        set(fh,'currentaxes',ax(1))
        shadedErrorBar(cin_t(good_cint),cin_hit_psth,cin_hit_psth_sem,...
            {'linestyle','-','color',this_color,'linewidth',1,'parent', ax(1)},1);
        set(fh,'currentaxes',ax(2))
        shadedErrorBar(cout_t(good_coutt),cout_hit_psth, cout_hit_psth_sem,...
            {'linestyle','-','color',this_color,'linewidth',1,'parent', ax(2)},1);
    else
        plot(ax(1),cin_t(good_cint),cin_hit_psth,'color',this_color,'linewidth',1.5);
        plot(ax(2),cout_t(good_coutt),cout_hit_psth,'color',this_color,'linewidth',1.5);
    end
    
end
linkaxes(ax,'y')

ylim(ax(1),([floor(min(psths(:))*10)/10 ceil(max(psths(:))*10)/10]))

xlim(ax(1),cintrange)
xlim(ax(2),couttrange)

%ax(2).YColor = 'w';
ax(1).TickDir = 'out';
ax(2).TickDir = 'out';

ylabel(ax(1), 'firing rate (spks/s)')
xlabel(ax(1), 'time from stim onset (s)')
xlabel(ax(2), 'from movement (s)')

%%
colormap(cm);
% cb = colorbar('north');
% cb.Position = cb.Position + [.025 .1 -.15 -.05];
cb = [];
cb = colorbar('east');
cb.Position = cb.Position + [.1 .15 -.0 -.3];
ax(2).YColor = 'w';
title(cb,cblab);
set(cb,'ytick',[]);
% drawnow
% ax1pos  = get(ax(1),'position');
% ax2pos  = get(ax(2),'position');
% posrat  = diff(cintrange) ./ diff(couttrange); 
% set(ax(2),'position',[ax2pos(1) ax1pos(2) ax1pos(3)/posrat ax2pos(4)])

%set(cb,'xtick',0:1,'xticklabel',b([1:2:end])) 
%%


if ncells == 1
    group_name = num2str(cell_to_plot);
else
    group_name = 'group';
end
print(fh, fullfile(dp.psth_fig_dir, ['cell_' group_name ]),...
    '-dsvg', '-painters')

