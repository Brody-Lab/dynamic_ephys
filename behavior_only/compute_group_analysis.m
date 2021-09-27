function [analysis] = compute_group_analysis(F);

%NLL_to_L = @(NLL) exp(-NLL);
compute_AIC = @(NLL, num_params, num_trials) 2*num_params + 2*NLL;
compute_BIC = @(NLL, num_params, num_trials) log(num_trials)*num_params + 2*NLL;

analysis = zeros(length(F), 2);
for i=1:length(F);
    if isstruct(F{i})
       analysis(i,1) = F{i}.f;
       analysis(i,2) = F{i}.fpt;       
       analysis(i,3) = compute_AIC(F{i}.f,  length(F{i}.final), F{i}.nt);
       analysis(i,4) = compute_BIC(F{i}.f,  length(F{i}.final), F{i}.nt);
    else
        error('rat doesnt have analysis')
    end
end

%%%%%%
if 0
%C = zeros(3,3);
VARS = zeros(1,8);
for i=1:length(F)
    if isstruct(F{i})
        F{i}.C = inv(F{i}.auto_hessian);
        vars = F{i}.C(8,:)./F{i}.C(8,8)
        VARS = VARS + vars;
        keyboard
    end
end


VARS = VARS./length(F);

end
