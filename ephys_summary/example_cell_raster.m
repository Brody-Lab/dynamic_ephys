function [fh, ax] = example_cell_raster(varargin)
p = inputParser;
addParameter(p,'cells',[ 16857 17784 18181 ]);
addParameter(p,'meta',0)
addParameter(p,'flip',0)
addParameter(p,'pause',0)
addParameter(p,'top_color','')
addParameter(p,'bot_color','')
addParameter(p,'max_ntrials',Inf)
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
meta = p.Results.meta;
flip = p.Results.flip;
fig_num     = p.Results.fig_num;
bin_size    = p.Results.bin_size;
max_nt      = p.Results.max_ntrials;

dp = set_dyn_path;
%%
repack = p.Results.repack;

min_t = p.Results.min_t;

[d, array_data, vec_data]  = dyn_cell_packager(cellid,'repack',repack);

j   = find(ismember(d.align_strs,coutstr));
i   = find(ismember(d.align_strs,cinstr));
k   = strmatch('stimend',d.align_strs,'exact');
if meta & p.Results.flip == 0
    flip = ~d.prefsideright{k};
end

T           = d.trials.T; %#ok<*AGROW>
sides  = 2.*d.trials.genEndState-1;
poke_r = d.trials.rat_dir==1;
gamma  = d.trials.gamma;

if flip==1
    gamma  = -gamma;
    sides  = -this_sides;
    poke_r = ~this_poke_r;
end

hit     = d.trials.hit==1;
end_state_s = d.trials.end_state_dur(:);

nt = length(array_data);

[cout_fr, cout_t] = make_rate_functions(cellid, ...
    'array_data', array_data, 'vec_data', vec_data,...
    'bin_size', bin_size, d.align_args{j}{:},...
    'krn_type', 'raster');

[cin_fr, cin_t] = make_rate_functions(cellid, ...
    'array_data', array_data, 'vec_data', vec_data,...
    'bin_size', bin_size, d.align_args{i}{:},...
    'krn_type', 'raster');
%%
good = T > min_t;

if sum(good) > max_nt
    [~, subsamprank] = sort(rand(size(good)));
    subsampind  = find(subsamprank <= max_nt);
    subsamp     = false(size(good));
    subsamp(subsampind) = true;
    good = good & subsamp;   
end
nt = sum(good);
go_r    = d.trials.rat_dir==1;
rng(0);
[go_r_sort, sort_ind] = sort(1000*go_r+T);
tmarker = vec_data.stim_dur;
go_r_sort   = go_r_sort(good); 
sort_ind    = sort_ind(good);
cin_raster  = cin_fr(sort_ind,:);
cin_raster(isnan(cin_raster)) = 0;
[i j]   = find(cin_raster);
ind_r   = go_r_sort(i)>=1000;
ind_l   = go_r_sort(i)<1000;

trial_markers_ = sum(ind_l);

fh = figure(fig_num); clf

ax(1) = subplot(121);hold(ax(1),'on');
ax(2) = subplot(122);hold(ax(2),'on');
axis(ax(1),'ij')
axis(ax(2),'ij')
xlim(ax(1),cintrange)
xlim(ax(2),couttrange)

plot(ax(1), cin_t(j(ind_r)), i(ind_r),'.',...
    'color',dp.right_color,'markersize',1)
plot(ax(1), cin_t(j(ind_l)), i(ind_l),'.',...
    'color',dp.left_color,'markersize',.1)
plot(ax(1), tmarker(sort_ind), 1:nt,'.r',...
    'markersize', 6)

%%
cout_raster = cout_fr(sort_ind,:);
cout_raster(isnan(cout_raster)) = 0;
[i j]   = find(cout_raster);
ind_r   = go_r_sort(i)>=1000;
ind_l   = go_r_sort(i)<1000;

plot(ax(2), cout_t(j(ind_r)),i(ind_r),'.',...
    'color',dp.right_color,'markersize',.1)

plot(ax(2), cout_t(j(ind_l)),i(ind_l),'.',...
    'color',dp.left_color,'markersize',.1)

tmarker = vec_data.stim_start-vec_data.cpoke_out;
plot(ax(2), tmarker(sort_ind), 1:nt,'.r',...
    'markersize',6)

axis ij

midl      = find(diff(go_r_sort<1000))+.5;

plot(ax(1), cintrange, midl*[1 1], 'k')
plot(ax(2), couttrange, midl*[1 1], 'k')

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

