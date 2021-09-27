function [mean_fr, var_fr, base_fr, dprime, last_mean_fr, last_var_fr, last_dprime] = compute_average_firing_rate(array_data,vec_data)

% get spike rates over entire trial
fr      = zeros(1, length(array_data));
basefr  = zeros(1, length(array_data));
endfr   = zeros(1, length(array_data));
for i=1:length(array_data)
    trial_spikes = array_data(i).spikes >= 0  & array_data(i).spikes <= array_data(i).stim_end;
    non_trial_spikes = array_data(i).spikes <= 0 & array_data(i).spikes >= -2;
    end_spikes = array_data(i).spikes >= (array_data(i).stim_end-.2)  & array_data(i).spikes <= array_data(i).stim_end;

    T       = array_data(i).stim_end;
    fr(i)   = sum(trial_spikes)./T;
    basefr(i) = sum(non_trial_spikes)./2;
    endfr(i)= sum(end_spikes)./0.2;
end

% compute some basic statistics
mean_fr = nanmean(fr);
var_fr  = var(fr);
R_fr    = fr(vec_data.pokedR);
L_fr    = fr(~vec_data.pokedR);
dprime  = (mean(R_fr) - mean(L_fr))./sqrt(.5*(var(R_fr)+var(L_fr)));
base_fr = nanmean(basefr);

last_mean_fr = nanmean(endfr);
last_var_fr  = var(endfr);
last_R_fr    = endfr(vec_data.pokedR);
last_L_fr    = endfr(~vec_data.pokedR);
last_dprime  = (mean(last_R_fr) - mean(last_L_fr))./sqrt(.5*(var(last_R_fr)+var(last_L_fr)));

scramble = randperm(length(vec_data.pokedR));
Rs = fr(scramble(vec_data.pokedR));
Ls = fr(scramble(~vec_data.pokedR));
s_dprime  = (mean(Rs) - mean(Ls))./sqrt(.5*(var(Rs)+var(Ls)));

scramble = randperm(length(vec_data.pokedR));
Rs = endfr(scramble(vec_data.pokedR));
Ls = endfr(scramble(~vec_data.pokedR));
s_last_dprime  = (mean(Rs) - mean(Ls))./sqrt(.5*(var(Rs)+var(Ls)));


