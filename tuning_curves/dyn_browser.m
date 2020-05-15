function [] = dyn_browser(cellid)

if nargin < 1 | isempty(cellid)
    % get list of cells in the MGO experiment
[~, cellid] = dyn_get_cells();
end

opts.bin_size=0.01;
opts.kx=-1:opts.bin_size:1;
opts.krn=normpdf(opts.kx,0, .1);  % kernel is a half-gaussian with sd of 100 ms.
opts.krn(opts.kx<0)=0;
opts.pre=2.5;
opts.post=2;
opts.wndw=1;
opts.xax=-opts.pre:opts.bin_size:(opts.post-opts.bin_size);
S.opts = opts;

% load the session data for those cells 
S.CD=get_celldata(cellid);
S.cellid=S.CD.cellid;
S.SD=[];
S.SD=get_sessdata(S.CD.sessid);

S.ncells = length(S.cellid);

S.fh = figure('units','pixels',...
              'position',[400 400 1800 900],...
              'menubar','none',...
              'name','slider_plot',...
              'numbertitle','off',...
              'resize','off');    

S.cin_psth_ax = axes('unit','pix',...
            'position',[50 100 500 350]);
        
S.cin_rast_ax = axes('unit','pix',...
            'position',[50 500 500 350]);
        
S.switch_psth_ax = axes('unit','pix',...
    'position',[600 100 500 350]);

S.switch_rast_ax = axes('unit','pix',...
    'position',[600 500 500 350]);

S.cout_psth_ax = axes('unit','pix',...
    'position',[1150 100 500 350]);

S.cout_rast_ax = axes('unit','pix',...
    'position',[1150 500 500 350]);

%linkaxes([S.cin_psth_ax,S.switch_psth_ax, S.cout_psth_ax],'y')
     
S.sl = uicontrol('style','slide',...
                 'unit','pix',...
                 'position',[50 10 800 30],...
                 'min',1,'max',length(S.cellid),'val',1,...
                 'sliderstep',[1/length(S.cellid) 1/length(S.cellid)],...
                 'callback',{@sl_call,S});  
             
function [] = sl_call(varargin)
    % Callback for the slider.
    [h,S] = varargin{[1,3]};  % calling handle and data structure.
    %set(S.LN,'ydata',S.x.^get(h,'value'))
    %S.LN = plot(S.x,S.x,'r');
    % update slider
    cla(S.cout_rast_ax)
    cla(S.cout_psth_ax)
    cla(S.cin_rast_ax)
    cla(S.cin_psth_ax)
    cla(S.switch_rast_ax)
    cla(S.switch_psth_ax)
    Value = round(get(h, 'Value'));
    set(h, 'Value', Value);
    disp(get(h,'value'))
    
    %%
    
    this_cell = S.cellid(Value);
    d = dyn_cell_packager(this_cell);
    sessid = bdata('select sessid from cells where cellid={S}',S.cellid(Value));
    SD = get_sessdata(sessid);
    bd =  [SD.pd{1}.bupsdata{:}];
    
    has_switch = ~cellfun(@isempty, d.trials.genSwitchTimes);
    last_switch = nan(size(has_switch));
    last_switch(has_switch) = cellfun(@(x) x(end), d.trials.genSwitchTimes(has_switch));
    d.trials.last_switch = last_switch;
    
    d.trials.bup_diff;
    d.trials.correct_dir;
    d.trials.rat_dir;
    d.trials.hit;
    
    stim_dur = d.trials.cpoke_end - d.trials.stim_start;
    cout_from_in = d.trials.cpoke_out - d.trials.stim_start;
    cin_from_cout = d.trials.stim_start - d.trials.cpoke_out;
    cout_from_switch = d.trials.cpoke_out - d.trials.last_switch;
    cin_from_switch = d.trials.stim_start - d.trials.last_switch;
    [stim_dur_sorted, stim_dur_sort_order] = sort(stim_dur);
    [cout_sorted, cout_sort_order] = sort(cout_from_in);
    
    stim_rast = d.frate{9};
    stim_rast_t = d.frate_t{9};
    cout_rast = d.frate{8};
    cout_rast_t = d.frate_t{8};
    switch_rast = d.frate{15};
    switch_rast_t = d.frate_t{15};
    
    choices = d.trials.rat_dir;
    [choice_sorted, choice_sort_order] = sort(choices);
    
    rast_sort_order = flipud(stim_dur_sort_order);
    rast_sort_order = flipud(cout_sort_order);
    rast_sort_order = rast_sort_order(choice_sort_order);
    
    easyR = d.trials.bup_diff > 30;
    easyL = d.trials.bup_diff < -30;
    hardR = d.trials.bup_diff < 30  & d.trials.bup_diff > 0;
    hardL = d.trials.bup_diff > -30 & d.trials.bup_diff < 0;
    co = colormapRedBlue(2).^.8.*.75;
    co(3,:) = [];
    co = flipud(co);
    

    % cin/stim raster
    imagesc(S.cin_rast_ax, (stim_rast(rast_sort_order,:)), 'x', stim_rast_t,[0 5]);
    hold(S.cin_rast_ax,'on');
    plot(S.cin_rast_ax,stim_dur(rast_sort_order),1:length(rast_sort_order),'k.')
    plot(S.cin_rast_ax,cout_from_in(rast_sort_order),...
        1:length(rast_sort_order),'b.')
    xlabel(S.cin_rast_ax, 'time from stim on/nose in (s)')
    colormap(S.cin_rast_ax,(flipud([colormapReds])))
    title(S.cin_rast_ax, [num2str(this_cell) ' side = ' num2str(d.brainsideright)])
    
    % cout raster
    imagesc(S.cout_rast_ax, (cout_rast(rast_sort_order,:)), 'x', cout_rast_t,[0 5]);
    hold(S.cout_rast_ax,'on');
    plot(S.cout_rast_ax,cin_from_cout(rast_sort_order),1:length(rast_sort_order),'k.')
    xlabel(S.cout_rast_ax, 'time from center out (s)')
    colormap(S.cout_rast_ax,(flipud([colormapReds])))
    
    % switch raster
    imagesc(S.switch_rast_ax, switch_rast(rast_sort_order,:), 'x', switch_rast_t,[0 5]);
    hold(S.switch_rast_ax,'on');
    plot(S.switch_rast_ax,cin_from_switch(rast_sort_order),1:length(rast_sort_order),'k.')
    plot(S.switch_rast_ax,cout_from_switch(rast_sort_order),...
        1:length(rast_sort_order),'bo')
    xlabel(S.switch_rast_ax, 'time from last switch (s)')
    colormap(S.switch_rast_ax,(flipud([colormapReds])))
    
    % cin/stim psth
    plot(S.cin_psth_ax,stim_rast_t, nanmean(stim_rast(easyR,:)),'linewidth',2,'color',co(1,:));
    hold(S.cin_psth_ax, 'on');
    plot(S.cin_psth_ax,stim_rast_t, nanmean(stim_rast(hardR,:)),'linewidth',2,'color',co(2,:));
    plot(S.cin_psth_ax,stim_rast_t, nanmean(stim_rast(hardL,:)),'linewidth',2,'color',co(3,:));
    plot(S.cin_psth_ax,stim_rast_t, nanmean(stim_rast(easyL,:)),'linewidth',2,'color',co(4,:));
    plot(S.cin_psth_ax,[0 0], get(S.cin_psth_ax,'ylim'),'k')
    xlabel(S.cin_psth_ax, 'time from stim on/nose in (s)')
    ylabel(S.cin_psth_ax,'spks/s')
    set(S.cin_psth_ax,'fontsize',17)
    
    % switch psth
    
    demean_switch_rast_endR = switch_rast(d.trials.genEndState==1&d.trials.hit,:) - nanmean(switch_rast);
        
    demn_switch_rast_endL = switch_rast(d.trials.genEndState==0&d.trials.hit,:) - nanmean(switch_rast);
    plot(S.switch_psth_ax,switch_rast_t, nanmean(demean_switch_rast_endR),'linewidth',2,'color',co(1,:));
    hold(S.switch_psth_ax, 'on');
    plot(S.switch_psth_ax,switch_rast_t, nanmean(demean_switch_rast_endR)-nansem(demean_switch_rast_endR),':','linewidth',2,'color',co(1,:));
    plot(S.switch_psth_ax,switch_rast_t, nanmean(demean_switch_rast_endR)+nansem(demean_switch_rast_endR),':','linewidth',2,'color',co(1,:));
    plot(S.switch_psth_ax,switch_rast_t, nanmean(demn_switch_rast_endL),'linewidth',2,'color',co(4,:));
    plot(S.switch_psth_ax,switch_rast_t, nanmean(demn_switch_rast_endL)-nansem(demn_switch_rast_endL),':','linewidth',2,'color',co(4,:));
    plot(S.switch_psth_ax,switch_rast_t, nanmean(demn_switch_rast_endL)+nansem(demn_switch_rast_endL),':','linewidth',2,'color',co(4,:));
    plot(S.switch_psth_ax,[0 0], get(S.switch_psth_ax,'ylim'),'k')
    xlabel(S.switch_psth_ax, 'time from last switch (s)')
    ylabel(S.switch_psth_ax,'spks/s')
    set(S.switch_psth_ax,'fontsize',17)
    
    % cout psth
    plot(S.cout_psth_ax,cout_rast_t, nanmean(cout_rast(easyR,:)),'linewidth',2,'color',co(1,:));
    hold(S.cout_psth_ax, 'on');
    plot(S.cout_psth_ax,cout_rast_t, nanmean(cout_rast(hardR,:)),'linewidth',2,'color',co(2,:));
    plot(S.cout_psth_ax,cout_rast_t, nanmean(cout_rast(hardL,:)),'linewidth',2,'color',co(3,:));
    plot(S.cout_psth_ax,cout_rast_t, nanmean(cout_rast(easyL,:)),'linewidth',2,'color',co(4,:));
    plot(S.cout_psth_ax,[0 0], get(S.cout_psth_ax,'ylim'),'k')
    xlabel(S.cout_psth_ax, 'time from center out (s)')
    ylabel(S.cout_psth_ax,'spks/s')
    set(S.cout_psth_ax,'fontsize',17)
    

    
%%
% should make a button to show the waveform