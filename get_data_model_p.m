function [model, constant_x, allsame] = get_data_model_p(data, vec_data)

dp = set_dyn_path;
model_p_dir = dp.model_mean_dir;

model_posterior_fn = fullfile(model_p_dir, ...
    sprintf('model_posterior_%i.mat',data.sessid));

% load model posterior for all trials
if ~exist(model_posterior_fn,'file')
    model = compute_model_mean(data.ratname, data.sessid);
else
    load(model_posterior_fn, 'model');
end
% Need to filter trials from model to match data variable
model   = model(vec_data.good);
ntrials = length(vec_data.good);
% Check for variable a-binning, which happens with large sensory noise to save space
numabins = length(model(1).posterior.avals);
constant_x = model(1).posterior.avals;
allsame = true;
for ti=2:ntrials
    if length(model(ti).posterior.avals) ~= numabins
        allsame=false;
        if length(model(ti).posterior.avals) > numabins
            numabins=length(model(ti).posterior.avals);
            constant_x = -10:0.1:10;
        end
    end
end