function [rate_one, rate_zero, nt, choice, sides, varargout] = build_rates(S,p)

% Check how many trials we have
if sum(S.n_done(S.haz == p.haz)) > p.nt_max
    nt = p.nt_max;
else
    nt = sum(S.n_done(S.haz == p.haz));
end


% Set up variables
sides       = NaN(nt,1);        % Vector of correct sides
noise_sides = NaN(nt,1);        % Vector of noisy agent correct sides
model_sides = NaN(nt,1);        % Vector of model responses
hits        = NaN(nt,1);        % Vector for rat hits
rate_one    = NaN(p.excess.steps,nt);    % Click rates in state one
rate_zero   = NaN(p.excess.steps,nt);    % Click rates in state zero
nt_done     = 0;                % Track what trial number we are on
trial_used  = false(nt,1);
trialids    = nan(nt,1);
modelGoR = compute_model_behavior(S.ratname{1}, p);
assert(length(modelGoR) == nt)

%%

% Organize Data
j = length(S.sessid);
while j > 0 & nt_done <= nt% see how many sessions we have left
    if S.haz(j) == p.haz   % Skip sessions without the correct hazard rate
        sess_trial_start = 1+sum(S.n_done(1:j))-S.n_done(j);
        for k=1:S.n_done(j)   % Iterate over trials
            trial_ind = sess_trial_start+k-1;
            
            % Check for number of trials done
            if nt_done >= nt;
                break
            end
            if S.pd{j}.bupsdata{k}.T >= p.excess.minT     % Ensure trial duration was long enough
                nt_done = nt_done +1;               % Keep track of which trial number we are on
                trialids(nt_done) = trial_ind;
                if trial_used(trial_ind)==1
                    keyboard
                end
                trial_used(trial_ind) = 1;
                lc = zeros(p.excess.full_steps, 1);               % Left Clicks
                rc = zeros(p.excess.full_steps, 1);               % Right Clicks
                %%%% DOWN SAMPLING CHANGE
                %            lc = zeros(p.excess.steps, 1);               % Left Clicks
                %            rc = zeros(p.excess.steps, 1);               % Right Clicks
                
                T = S.pd{j}.bupsdata{k}.T;               % Duration of this trial
                left = S.pd{j}.bupsdata{k}.left;         % Grab left clicks
                left(find(left - T + p.excess.minT - p.excess.dt < 0)) = [];  % Exclude clicks before time frame
                left = left - (T - p.excess.minT);                   % shift click times to minimum time frame
                left = round(left.*(1/p.excess.dt));                 % round, even though will be integers
                lc(left) = 1;                               % Set times of left clicks
                
                right = S.pd{j}.bupsdata{k}.right;       % Grab right clicks
                right(find(right - T + p.excess.minT - p.excess.dt < 0)) = [];% Exclude clicks before time frame
                right = right - (T - p.excess.minT);                 % shift click times to minimum time frame
                right = round(right.*(1/p.excess.dt));               % round, even though will be integers
                rc(right) = 1;                              % Set times of right clicks
                
                % compute click rates
                temp = conv(double(rc(:)),p.excess.window) - conv(double(lc(:)),p.excess.window);
                %%%% DOWN SAMPLING CHANGE
                temp = temp(1:p.excess.full_steps);
                bindex = ceil(((1:p.excess.full_steps)/p.excess.ds));
                tempds = accumarray(bindex',temp');
                %            tempds = tempds./p.excess.ds; average vs sum doesn't matter
                rate_one(:, nt_done) = tempds;
                rate_zero(:,nt_done) = tempds;
                %%            rate_one(:, nt_done) = temp(1:p.excess.ds:p.excess.full_steps);
                %%            rate_zero(:,nt_done) = temp(1:p.excess.ds:p.excess.full_steps);
                %            rate_one(:, nt_done) = temp(1:p.excess.steps);
                %            rate_zero(:,nt_done) = temp(1:p.excess.steps);
                
                % compute generative state
                gs = compute_state(T, S.pd{j}.bupsdata{k}.genSwitchTimes, S.pd{j}.bupsdata{k}.genEndState,p.excess.dt);
                %%%% DOWN SAMPLING CHANGE
                %           gs = gs(round((T-p.excess.minT)*(1/p.excess.dt))+1:end);
                gs = gs(end-p.excess.full_steps+1:end); % shorten to frame
                gs = gs(1:p.excess.ds:p.excess.full_steps); % downsample
                
                % Click rate is NaN when we are in the other state
                rate_one(~logical(gs),nt_done) = NaN;     % NaN out generative state 0 time points
                rate_zero(logical(gs),nt_done) = NaN;     % NaN out generative state 1 time points
                
                % Save track of hits and sides
                sides(nt_done) = S.pd{j}.sides(k);
                hits(nt_done) = S.pd{j}.hits(k);
                if p.optimal.compute
                    noise_sides(nt_done) = S.pd{j}.bupsdata{k}.noiseAnswer;
                end
                if p.compute_model
                    model_sides(nt_done) = modelGoR(trial_ind);
                end
                
            end     % End trial duration check
        end     % End iterate over trials
    end         % End skip sessions with wrong hazard rate
    j = j-1;    % Increase the iteration over sessions, regardless of whether hazard matches
end             % End iterate over sessions


% Find Median click rate for each of the generative states
if p.excess.mean == 1 | p.excess.mean == 2
    m_one  = nanmean(rate_one(:,:),2);
    m_zero = nanmean(rate_zero(:,:),2);
elseif p.excess.mean == 0
    error('WARNING, WHY ARE YOU USING MEDIAN???? MAJOR ANALYSIS PROBLEMS')
    m_one  = nanmedian(rate_one(:,:),2);
    m_zero = nanmedian(rate_zero(:,:),2);
else
    error('Bad Parameter value for p.excess.mean: options 0,1,2')
end

% Normalize local click rates
for k=1:nt_done
    rate_one(:,k) = rate_one(:,k) - m_one;
    rate_zero(:,k) = rate_zero(:,k) - m_zero;
end

% Clean up sides and hits with respect to nt_done
sides = sides(1:nt_done);
hits = hits(1:nt_done);
if p.optimal.compute
    noise_sides = logical(noise_sides(1:nt_done));
end
if p.compute_model
    model_sides = logical(model_sides(1:nt_done));
end

% Compute Rat's Choice
sides  = sides =='r';
choice = (sides & hits) | (~sides & ~hits); % 1= right, 0 = left


% add 5% lapse rate
if 0
    dex = rand(nt_done,1) < .10;
    choice(dex) = ~choice(dex);
end

if p.optimal.compute & p.compute_model
    varargout{1} = noise_sides;
    varargout{2} = model_sides;
elseif p.optimal.compute & ~p.compute_model
    varargout{1} = noise_sides;
    varargout{2} = [];
elseif ~p.optimal.compute & p.compute_model
    varargout{1} = [];
    varargout{2} = model_sides;
end

