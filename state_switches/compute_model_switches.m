function [array_data, model_mean] = compute_model_switches(array_data, sessid, varargin)
p = inputParser;
addParameter(p,'model_dir',[])
addParameter(p,'model_smooth_wdw',[])
addParameter(p,'remove_initial_choice',1)
addParameter(p,'eval_dt',1e-3);
addParameter(p,'strength_window',.1);
addParameter(p,'clear_bad_strengths',1);
addParameter(p,'bad_strength',0);
addParameter(p,'fit_line',1);
addParameter(p,'change_bounds',[0 0]);
addParameter(p,'t_buffers',[0 0]);

parse(p,varargin{:});
p = p.Results;

dp = set_dyn_path;
if isempty(p.model_dir)
    model_dir = dp.model_mean_dir;
else
    model_dir = p.model_dir;
end
assert(diff(p.change_bounds) >= 0)
ratname = bdata('select ratname from sessions where sessid={S}',sessid);
ratname = ratname{1};

model_mean_fn   = sprintf('model_mean_%i.mat',sessid);
model_mean_file = fullfile(model_dir, model_mean_fn);
% load the model mean trajectory

which_trials = [array_data.trialnum];

if exist(model_mean_file,'file')
    m = load(model_mean_file);
    model_mean = m.model_mean(which_trials);
else
    compute_model_mean(ratname, sessid)
    m = load(model_mean_file);
    model_mean = m.model_mean(which_trials);
end
% smooth the mean model trajectory
if ~isempty(p.model_smooth_wdw)
    for j = 1:length(model_mean)
        model_mean(j).raw_mean = model_mean(j).mean;
        model_mean(j).mean = movmean(model_mean(j).mean, p.model_smooth_wdw);
    end
end
% compute the model switch times
array_data = compute_model_state(array_data, model_mean);

array_data = quantify_model_state_strength(array_data,'eval_dt',p.eval_dt,...
    'strength_window', p.strength_window, 'fit_line', p.fit_line);

array_data = smooth_model_state(array_data,...
    'remove_initial_choice', p.remove_initial_choice,...
    'clear_bad_strengths',p.clear_bad_strengths,'bad_strength',p.bad_strength,...
    'change_bounds',p.change_bounds, 't_buffers',p.t_buffers);

array_data = get_model_state_durs(array_data);
end