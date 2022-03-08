function [S] = load_data(rat,p)

% log status
x = ['Loading data for ' rat];
disp(x);

% Pull session data from server
[sessids, n_done] = bdata(['select sessid, n_done_trials from sessions '...
 'where ratname regexp "{S}" and protocol="PBupsDyn"'],char(rat) );

if length(sessids) > p.ns
    sessids = sessids(end-p.ns:end);
    n_done  = n_done(end-p.ns:end);
end

% Pull data 
S = get_sessdata(sessids);
if ~isempty(sessids)
    for j=1:length(S.sessid)
        if sum(isnan(S.pd{j}.hits))>0
            S.pd{j}.hits(isnan(S.pd{j}.hits)) = 0;
        end
        S.acc(j) =  sum(S.pd{j}.hits) / length(S.pd{j}.hits);
        S.vio(j) =  sum(S.pd{j}.violations) / length(S.pd{j}.violations);
        S.haz(j) =  S.pd{j}.bupsdata{end}.hazard;
        S.r1(j)  =  S.pd{j}.bupsdata{end}.rate_1;
        S.r2(j)  =  S.pd{j}.bupsdata{end}.rate_2;
        S.pd{j}.stimdata = []; % Clear stim data
        S.n_done(j) = n_done(j);

    end
    S.peh = [];
end


% Remove all violations (including c-poke)
if p.clearViolations
    for i=1:length(S.sessid)
        % need to document change in analysis...excess click rates changing????
        % need to separate out different types of violations?
        dex =  logical(S.pd{i}.violations);
        S.pd{i}.hits(dex)       = [];
        S.pd{i}.violations(dex) = [];
        S.pd{i}.sounds(dex)     = [];
        S.pd{i}.samples(dex)    = [];
        S.pd{i}.memory_gaps(dex)= [];
        S.pd{i}.sides(dex)      = [];
        S.pd{i}.rt_task(dex)    = [];
        S.pd{i}.stims(dex)      = [];
        S.pd{i}.bupsdata        = {S.pd{i}.bupsdata{~dex}}';
        S.pd{i}.n_left(dex)     = [];
        S.pd{i}.n_right(dex)    = [];
    end
else
    % Just remove extra trial from bcontrol preparing next trial
    for i=1:length(S.sessid)
        S.pd{i}.bupsdata = {S.pd{i}.bupsdata{1:end-1}}';
    end
end


% compute behavior of a linear agent
if p.optimal.compute
    disp('Computing behavior of optimal linear agent (takes a while)')
    for i=1:length(sessids)
        if length(p.haz_compare) > 1
            dex = find(p.haz_compare == S.haz(i));
            if isempty(dex)
                lambda = p.optimal.lambda;
            else
                lambda = p.optimal.lambda_compare(dex);
            end
        else
            lambda = p.optimal.lambda;
        end
        for j=1:length(S.pd{i}.hits)
            S.pd{i}.bupsdata{j}.noiseAnswer = compute_linear_agent(S.pd{i}.bupsdata{j},p,lambda);
        end
    end
end

for i=1:length(S.sessid)
    if ~isempty(S.pd{i}.bupsdata)
    if ~isfield(S.pd{i}.bupsdata{end},'hazardBarrier')
        for j=1:length(S.pd{i}.hits)
            S.pd{i}.bupsdata{j}.hazardBarrier = 0;
        end
    end
    end
end


% % Load model fits and behavior
% if p.model.compute
% %    error('Model comparison not implemented yet!')
%     % See if model has been fit to this rat
%     filename = [p.datapath_root p.modelfit '/fit_' rat '.mat'];
%     if exist(filename) == 2
%         % Load model parameters       
%         disp(['Loading model parameters for rat: ' rat]);
%         load(filename,'history');
%         p.model.params = history.x(end,:);
% 
%         %load model responses that we already have computed. 
%         %a struct with sessid, and a vector of model R/L   
%         filename = [p.datapath_root rat p.model_data '.mat'];
%         if exist(filename) == 2
%             disp('Loading existing model responses')
%             load(filename); 
%         else
%             disp('No model responses found')
%             modelB = struct('sessid',[],'goR',[]);
%             modelB(1) = [];
%         end  
% 
%         modelSess = [modelB.sessid];      
%         update = 0;
%         disp('Checking for sessions that need model computations, this could take a while')
%         if length(modelSess) < p.ns
%         diff = p.ns - length(modelSess);
%         else
%         diff = 1;
%         end
%         disp(['Probably need to compute model behavior for ' num2str(diff) ' sessions']);
%         for i=1:length(S.pd) 
%             dex = find(modelSess == S.sessid(i));
%             if isempty(dex)
%             % If we don't have a session, compute model responses, and add to list
%                 disp(['Computing Model behavior for rat ' rat ' and session #' num2str(S.sessid(i))]);
%                 tic;
%                 update = update + 1;
%                 modelB(end+1).sessid = S.sessid(i);
%                 modelB(end).goR = compute_model_behavior([S.pd{i}.bupsdata],p);
%                 S.pd{i}.modelGoR = modelB(end).goR; 
%                 toc
%             elseif length(dex) == 1
%             % Iterate over sessids, incorporating model responses.
%             % check to make sure its the right length
%              %   disp(['Found model behavior for rat ' rat ' and session #' num2str(S.sessid(i))]); 
%                 if length(S.pd{i}.hits) == length(modelB(dex).goR)
%                     S.pd{i}.modelGoR = modelB(dex).goR;
%                 else
%                     error('Found model behavior that doesnt match rat data');
%                 end
%             else
%                 error('Sessid is found multiple times in model behavior database!')
%             end 
%         end
% 
%         % Save model responses if required 
%         if update > 0
%             filename = [p.datapath_root rat p.model_data '.mat'];
%             save(filename, 'modelB');
%             disp(['Computed Model responses for ' num2str(update) ' sessions']);
%         end
%     
%     else
%         disp(['No model fit found for rat: ' rat]);
%         % might want to put an empty field for S.pd{i}.modelGoR;
%     end
% end



% Check to see when the rat last trained, and report if its been more than 1 week
%detect_rat_failure(p,S);


