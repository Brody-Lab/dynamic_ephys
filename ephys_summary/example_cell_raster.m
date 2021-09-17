function [fh ax] = example_cell_psth(varargin)


p = inputParser;
addParameter(p,'cells',[ 16857 17784 18181 ]);
addParameter(p,'type','choice'); %evR, chrono
addParameter(p,'separate_hits',1);
addParameter(p,'meta',0)
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
addParameter(p, 'cinstr', 'stimstart-cout-mask');
addParameter(p, 'coutstr', 'cpokeout');
addParameter(p, 'do_print', 0);
addParameter(p, 'bin_size', .005);
parse(p,varargin{:});
cintrange = p.Results.cintrange;
couttrange = p.Results.couttrange;
cinstr = p.Results.cinstr;
coutstr = p.Results.coutstr;

cellid = p.Results.cells;
type = p.Results.type;
meta = p.Results.meta;
flip = p.Results.flip;
do_pause = p.Results.pause;
edges = p.Results.edges;
separate_hits = p.Results.separate_hits;
fig_num     = p.Results.fig_num;
bin_size    = p.Results.bin_size;


dp = set_dyn_path;
%%
repack = p.Results.repack;

[align_strs, align_args] = dyn_align_LUT(2);
i = find(strcmp(align_strs,cinstr));
j = find(strcmp(align_strs,coutstr));
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


[d, array_data, vec_data]  = dyn_cell_packager(cellid,'repack',repack);

end_tind = find(ismember(d.align_strs,coutstr));
start_tind  = find(ismember(d.align_strs,cinstr));
stimind = strmatch('stimend',d.align_strs,'exact');
if meta & p.Results.flip == 0
    flip = ~d.prefsideright{stimind};
end
T       = [T; d.trials.T]; %#ok<*AGROW>
this_sides = 2.*d.trials.genEndState-1;
this_poke_r = d.trials.rat_dir==1;
this_gamma  = d.trials.gamma;
if flip==1
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

nt = length(array_data);

[cout_fr, cout_t] = make_rate_functions(cellid, ...
    'array_data', array_data, 'vec_data', vec_data,...
    'bin_size', bin_size, align_args{j}{:},...
    'krn_type', 'raster');

[cin_fr, cin_t] = make_rate_functions(cellid, ...
    'array_data', array_data, 'vec_data', vec_data,...
    'bin_size', bin_size, align_args{i}{:},...
    'krn_type', 'raster');


%%

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
    %top = [.75 .75 .25];

    cm_bot = colormapLinear(bot, ceil(nbins/2));
    cm_top = colormapLinear(top, ceil(nbins/2));
    cm = [cm_bot(end:-1:2,:); cm_top(2:end,:)];
    
else
    cm = color_set(nbins);
end

err = hit == 0;

good = T > min_t;

psths = [];
 
good_cint   = cin_t > cintrange(1) & cin_t < cintrange(2);
good_coutt  = cout_t > couttrange(1) & cout_t < couttrange(2);

cin_raster = [];
cout_raster = [];

trial_markers = [];
%%
fh = figure(fig_num); clf
set(fh,'position',[15 5 6 3 ],'papersize',[6 3],'paperpositionmode','auto')

ax(1) = subplot(121);hold(ax(1),'on');
ax(2) = subplot(122);hold(ax(2),'on');
xlim(ax(1),cintrange)
xlim(ax(2),couttrange)
go_r    = d.trials.rat_dir==1;
go_l    = ~go_r;
[go_r_sort, sort_ind] = sort(go_r);

cin_raster = cin_fr(sort_ind,:);
cout_raster = cout_fr(sort_ind,:);
cin_raster(isnan(cin_raster)) = 0;
cout_raster(isnan(cout_raster)) = 0;
[i j]   = find(cin_raster);

ind_r   = go_r_sort(i);
ind_l   = ~go_r_sort(i);

nt      = length(go_r);

trial_markers_ = sum(ind_l);

tic
plot(ax(1), cin_t(j(ind_r)),i(ind_r),'.',...
    'color',dp.right_color,'markersize',.1)

plot(ax(1), cin_t(j(ind_l)),i(ind_l),'.',...
    'color',dp.left_color,'markersize',.1)
axis(ax(1),'ij')

[i j]   = find(cout_raster);
ind_r   = go_r_sort(i);
ind_l   = ~go_r_sort(i);
plot(ax(2), cout_t(j(ind_r)),i(ind_r),'.',...
    'color',dp.right_color,'markersize',.1)

plot(ax(2), cout_t(j(ind_l)),i(ind_l),'.',...
    'color',dp.left_color,'markersize',.1)
axis ij
toc

midl      = find(diff(go_r_sort))+.5;

plot(ax(1), cintrange, midl*[1 1], 'k')
plot(ax(2), couttrange, midl*[1 1], 'k')
%%
linkaxes(ax,'y')
colormap(flipud(bone))
%ylim(ax(1),([floor(min(psths(:))*10)/10 ceil(max(psths(:))*10)/10]))

set(ax(1), 'ytick', [100:100:1000])
set(ax(2), 'ytick', [100:100:1000])
ax(1).TickDir = 'out';
ax(2).TickDir = 'out';

ylabel(ax(1), 'Trial #')
xlabel(ax(1), 'Time from stim onset (s)')
xlabel(ax(2), 'from movement (s)')
plot(ax(1),[ 0 0], ylim,'k')
plot(ax(2),[ 0 0], ylim,'k')
ylim(ax(1),[0 nt+1])
ylim(ax(2),[0 nt+1])
%%

if p.Results.do_print
print(fh, fullfile(dp.psth_fig_dir, ['cell_' group_name '_raster']),...
    '-dsvg', '-painters')
end

