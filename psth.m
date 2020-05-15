function [x,y,yse,rx,ry,n] = psth(array_data,vec_data, time_field, trial_dex,p);

% set up smoothing filter
dt      = 5e-4;
fstd     = .1;
kernel  = exp(-.5.*(dt:dt:fstd*4).^2./fstd^2);
kernel  = kernel./sum(abs(kernel));

% filter trials
if ischar(trial_dex) & strcmp(trial_dex, 'all')
    trial_dex = logical(ones(1,length(vec_data.pokedR)));
end
 
% align spikes to reference events
[ref, ts]   = compute_ref(array_data,vec_data,time_field,trial_dex);

% convolve against each reference event
[y x r]   = spike_filter(ref, ts, kernel);

% Need to do anything special?
if p.sort_by_state == 0
    p.normalize = 0;
    reft = find(x==0);
    for i=1:length(array_data) 
        temp = ones(1,size(y,2));
        bad_state = array_data(i).gen_state(1:5:end) == 1;
        if strcmp(time_field,'stim_start')
            temp(reft+1:reft+length(bad_state)) = bad_state;
        else
            temp(reft+1-length(bad_state):reft) = bad_state;
        end
        y(i,logical(temp)) = NaN;
    end
elseif p.sort_by_state == 1
    reft = find(x==0);
    p.normalize = 0;
    for i=1:length(array_data)
        temp = ones(1,size(y,2));
        bad_state = array_data(i).gen_state(1:5:end) == 0;
        if strcmp(time_field, 'stim_start')
            temp(reft+1:reft+length(bad_state)) = bad_state;
        else
            temp(reft+1-length(bad_state):reft) = bad_state;   
        end
        y(i,logical(temp)) = NaN;
    end
end

if strcmp(p.sort_by_state_dur, 'all')
    p.normalize = 0;
    reft = find(x==0);
    durs_vec= 0:5e-4:2; 
    durs    = zeros(1,length(durs_vec));
    counts  = zeros(1,length(durs));
    for i=1:length(array_data)
        for j=1:5:length(array_data(i).gen_state_duration)
            ydex = reft + floor(j/5);
            durdex          = round(array_data(i).gen_state_duration(j)/5e-4);
            if durdex == 0; durdex =1; end;
            durs(durdex)    = durs(durdex) + y(i,ydex);
            counts(durdex)  = counts(durdex) + 1;
        end
    end
    durs = durs./counts;
    y = durs;
    x = durs_vec;
    yse = zeros(size(y));
    rx = [];
    ry = [];
    n = length(ref);
    return 
end

% compute mean and se
if p.normalize
    y       = y - y(:,find(x == 0));    
end

yse = nanstd(y)./sqrt(length(ref));
y = nanmean(y);

rx = [];
ry = [];
n = length(ref);

if p.raster
    num_spikes = sum(r(:));
    raster_x = nan(num_spikes*3,1);
    raster_y = nan(num_spikes*3,1);
    c = 1;
    for i=1:n
        these_spikes = x(find(r(i,:)));
        nc = c + length(these_spikes);
        sd = (c-1)*3+1;
        ed = (nc-1)*3;
        raster_x(sd:3:ed)        = these_spikes;
        raster_x(sd+1:3:ed+1)    = these_spikes;
        raster_x(sd+2:3:ed+2)   = nan;
        
        raster_y(sd:3:ed)        = i;
        raster_y(sd+1:3:ed+1)    = i+.8;
        raster_y(sd+2:3:ed+2)   = nan;
        c = nc;
    end
    rx = raster_x;
    ry = raster_y;
end    
