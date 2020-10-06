function [results] = compute_variance_explained(results, p,cellids, force)

    % If we weren't given a list of cells, compute the variance explained for all of them
    if nargin ==  1
        cellids = results.cellid; 
        force = 0;
        p.lag = 0;
    end
    if ~isfield(results, 'var_explained')
        results.var_explained = NaN(size(results.frm_time));
    end

    % Iterate through list of cells, and compute the variance explained
    for i=1:length(cellids)
        lag = p.lag.*10;
        if lag == 0;
            lagstr = '';
        else
            lagstr = ['_' num2str(lag)];
        end
        if exist(['/home/alex/Dropbox/spikes/data/var_explained/' num2str(cellids(i)) '_var_explained' lagstr '.mat'],'file') & ~force
            disp('variance explained already computed')
            load(['/home/alex/Dropbox/spikes/data/var_explained/' num2str(cellids(i)) '_var_explained' lagstr '.mat'])
            celldex = find(results.cellid == cellids(i));
            results.var_explained(:,celldex) = var_explained;
            results.var_explained_static(:,celldex) = var_explained_static;
        else
            disp('computing variance explained')
            [results, var_explained, var_explained_static] = compute_variance_explained_single(results, cellids(i),p);    
            save(['/home/alex/Dropbox/spikes/data/var_explained/' num2str(cellids(i)) '_var_explained' lagstr '.mat'],'var_explained', 'var_explained_static')
        end
    end
end


function [results, var_explained, var_explained_static] = compute_variance_explained_single(results, cellid,p)
    % For each trial I need 
    %   actual r(t)
    %   a(t)
    %   model predicted R = r(a,t)    
    
    % to get r(a,t), load from results matrix 
    celldex = find(results.cellid == cellid);
    if isempty(celldex)
        error('Couldnt find cell in cellid list')   
    end
    if isfield(results, 'tuning_cell')
        fga = results.tuning_cell(:,:,celldex) + nanmean(results.fga_cell(:,:,celldex),2) ; %%same as results.fga_aa_cell
        fga_mean = nanmean(results.fga_cell(:,:,celldex),2) ;
        if results.flipdex(celldex)
            fga = flipdim(fga,2);
        end

    else
        error('Couldnt find tuning curve, do you want to use the default fga_cell?')
    end
   
    % to get a(t), get cellid/sessid, then load model_mean
    data    = dyn_cell_packager(cellid);
    sessid  = data.sessid; 
    load(['/home/alex/Dropbox/spikes/model/model_mean_' num2str(sessid) '.mat'])
    % Need to filter trials from model to match data variable
    p.reload    = 0;
    p.ratname   = data.ratname;
    [~, vec_data, ~,~] = get_behavior_data('/home/alex/Dropbox/spikes/data/', data.cellid, data.sessid,p);
    model_mean = model_mean(vec_data.good);

    % to get r(t), I need to determine the alignment index
    alignment = 'stimstart-cout-mask';
    [align_strs, align_args] = dyn_align_LUT;
    align_ind = strmatch(alignment,align_strs,'exact');
    fr_all = data.frate{align_ind}; 
    ft = data.frate_t{align_ind};
    
    % compute model predicted R via r(a,t)
    R   = nan(size(fr_all));
    Rm  = nan(size(fr_all));
    t0s = results.t0s;
    if min(diff(t0s)) < min(diff(ft))
        error('Time bins for tuning curves are smaller than firing rate smoothed vector. Recompute firing rate.')
    end
    % Check to see if the time bins requested for the tuning curve extend beyond the firing rate vector
    if t0s(end) > ft(end)
        error('The last time point is beyond the firing rate time points. Decrease the time bins t0s, or recompute the firing rate')
    end
    
    % iterate through trials
    for i=1:size(R,1)
        fr = fr_all(i,:);
        t = model_mean(i).T;
        % iterate through timepoints
        for j=1:length(t0s);       
            t0 = t0s(j);
            t0_durs     = diff(t0s);
            if j == numel(t0s)
                this_dur    = 0;
            else
                this_dur    = t0_durs(j);
            end
            if j > 1
                last_dur = t0_durs(j - 1);
            else
                last_dur = 0;
            end
            if (this_dur/2 + t0) > ft(find(~isnan(fr),1,'last'))
                this_dur = 2*(ft(find(~isnan(fr),1,'last')) - t0);
            end
            all_time_indexes = t > (t0-last_dur/2) & t < (t0 + this_dur/2);
            mean_aval = mean(model_mean(i).mean(all_time_indexes));
            adex = bin_aval(results.dv_axis,mean_aval); %is dv axis centers or edges?
 
            % need to take mean firing rate within each time bin
            t0_durs     = diff(t0s);
            if j == numel(t0s)
                this_dur    = 0;
            else
                this_dur    = t0_durs(j);
            end
            if j > 1
                last_dur = t0_durs(j - 1);
            else
                last_dur = 0;
            end
            if (this_dur/2 + t0+lag) > ft(find(~isnan(fr),1,'last'))
                this_dur = 2*(ft(find(~isnan(fr),1,'last')) - t0-lag);
            end
            r0_mesh = interp1(ft, fr, [t0+lag-last_dur/2, t0+lag, t0+lag+this_dur/2]);
            r0      = mean(r0_mesh);% The mean firing rate during the time window
        
            % make sure you mask after stimulus firing
            include_fr = ~isnan(r0) & (t(end) -t0-lag  > 0);
            if include_fr
                R(i,j)  = r0; % store raw firing rate
                Rm(i,j) = fga(j,adex);
            end
        end
    end

    % compute variance explained
    Rmean = nanmean(R,1);
    total_var = nanmean((R-Rmean).^2,1);
    model_var = nanmean((R-Rm).^2,1);
    pVar = 1 - model_var./total_var;

    Rmean_static = nanmean(R(:));
    total_var_static = nanmean((R-Rmean_static).^2,1);
    model_var_static = nanmean((R-Rm).^2,1);
    pVar_static = 1 - model_var_static./total_var_static;

if 0
    figure(1);clf
    trial_res = 25; 
    subplot(3,3,1); hold on;
    plot(R(1:trial_res:end,:)')
    plot(Rmean, 'k','linewidth',2)
    xlabel('time')
    ylabel('data fr (hz)')
    x = ylim();
    subplot(3,3,2); hold on;
    plot(Rm(1:trial_res:end,:)')
    plot(Rmean, 'k','linewidth',2)
    plot(fga_mean, 'm','linewidth',2)
    ylim(x)
    xlabel('time')
    ylabel('model fr (hz)')
    subplot(3,3,5); hold on;
    plot(total_var,'k')
    plot(model_var,'r')
    xlabel('time')
    ylabel('Variance')
    subplot(3,3,4); hold on;
    plot((R(1:trial_res:end,:)-Rmean).^2')
    xlabel('time')
    ylabel('data fr. variance ')
    subplot(3,3,3); hold on;
    plot((R(1:trial_res:end,:)-Rm(1:trial_res:end,:)).^2')
    xlabel('time')
    ylabel('model fr. variance ')
    subplot(3,3,6); hold on;
    plot(pVar, 'k','linewidth',2)
    xlabel('time')
    ylabel('VE %')
    subplot(3,3,7); hold on;
    plot((R(1:trial_res:end,:)-Rmean_static).^2')
    xlabel('time')
    ylabel('data fr. variance ')
    subplot(3,3,8); hold on;
    plot(total_var_static,'b')
    plot(model_var_static,'r')
    xlabel('time')
    ylabel('Variance')
    subplot(3,3,9); hold on;
    plot(pVar_static, 'b','linewidth',2)
    xlabel('time')
    ylabel('VE %')
end
    results.var_explained(:,celldex) = pVar(1:length(t0s));
    var_explained = pVar(1:length(t0s));
    results.var_explained_static(:,celldex) = pVar_static(1:length(t0s));
    var_explained_static = pVar_static(1:length(t0s));
end




