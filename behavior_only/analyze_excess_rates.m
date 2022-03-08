function clicks = analyze_excess_rates(rat, S, p)

% log current step
disp(['Plotting Excess Click Rates for ' rat]);

% Analysis parameters
p.excess.win_dur    = 0.5;

if ~isfield(p.excess, 'OVERRIDE')
    p.excess.wStd       = 80;
    p.excess.ds         = 10;   
end 
%if p.haz == 0
%p.excess.wStd = 500;
%p.excess.ds = 5;
%end

p.excess.dt         = 0.0001;
p.excess.ww         = p.excess.wStd*4;
if p.haz == 0
    p.excess.mean       = 0;
else
    p.excess.mean       = 1;
end
p.excess.normalize  = 1;
p.excess.minT       = p.excess.win_dur + .05;
p.excess.window=exp(-.5.*(1:p.excess.ww).^2./p.excess.wStd^2)./sum(exp(-.5.*(1:p.excess.ww).^2./p.excess.wStd^2));
p.excess.full_steps = length(p.excess.dt:p.excess.dt:p.excess.minT);
p.excess.steps      = length(p.excess.dt:p.excess.dt*p.excess.ds:p.excess.minT);
p.compute_model     = 1;
% 
% % Compute behavior of linear model
% p.optimal.compute   = 0;
p.optimal.lambda    = -4.1;
S = compute_linear_agent_dataset(rat,S,p);

% Call helper function that compute click rates
warning('off','all')
clicks = compute_excess_rates(S, p);


% Save data for summary

save_dir = fullfile(p.data_dir, rat);
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
save_fn = fullfile(save_dir, ['excess' p.hazStr]);
save(save_fn, 'clicks');
disp(['Data saved to:  ' save_fn]);

% Turn warnings back on
warning('on','all')
