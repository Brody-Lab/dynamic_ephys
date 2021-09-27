function [model, constant_x, allsame] = get_data_model_p(data, which_trials, varargin)
% function [model, constant_x, allsame] = get_data_model_p(data, which_trials)
% returns the model posterior stuct for a given session for which_trials trials
% can take as input data either:
%   - filepath to model posterior
%   - sessid for finding model posterior file
%   - struct with fields sessid and ratname
%   - struct containing model
% model will be returned with only which_trials (or all trials if nargin == 1)
% output allsame tells if accumulator value is binned the same on each trial

if isstr(data)
    posterior_path = data;
elseif isnumeric(data)
    sessid  = data;
    ratname = bdata('select ratname from sessions where sessid={S}',sessid);
    ratname = ratname{1};
elseif isfield(data, 'sessid')
    sessid = data.sessid;
    ratname = data.ratname;
elseif isfield(data,'posterior')
    model = data;
end

if ~exist('model','var')
    if ~exist('posterior_path','var') 
        dp = set_dyn_path;
        model_p_dir = dp.model_mean_dir;
        posterior_path = fullfile(model_p_dir, ...
            sprintf('model_posterior_%i.mat',sessid));
    end
    % load model posterior for all trials
    if ~exist(posterior_path,'file')
        if ~exist('sessid','var')
            error('bad inputs');
        end
        model = compute_model_mean(ratname, sessid, varargin{:});
    else
        load(posterior_path, 'model');
    end
end

if nargin < 2
    which_trials = 1:length(model);
elseif isstruct(which_trials) % if it's a struct assume it's vec_data
    which_trials = which_trials.good;
end

% Need to filter trials from model to match data variable
model   = model(which_trials);
% Check for variable a-binning, which happens with large sensory noise to save space
constant_x  = model(1).posterior.avals;
navals      = cellfun(@(x) length(x.avals), ...
    {model.posterior});
allsame     = length(unique(navals)) == 1;

if ~allsame
    constant_x = -10:.1:10;
end
