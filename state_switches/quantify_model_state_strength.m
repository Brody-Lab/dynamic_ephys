function [array_data] = quantify_model_state_strength(array_data,varargin)
p = inputParser;
addParameter(p,'eval_dt',1e-3);
addParameter(p,'strength_window',.1);
addParameter(p,'fit_line',1);
parse(p,varargin{:});
p = p.Results;

% calculates the strength of each model state change.
% does so by calculating the mean model trajectory in the window...
% before and after each state change
nw = ceil(p.strength_window/p.eval_dt);
for i=1:length(array_data)
    for j=1:length(array_data(i).model_switch_to_0)
        dex = find(array_data(i).model_T == array_data(i).model_switch_to_0(j));    
        sd = dex - nw;
        ed = dex + nw;
        if sd < 1;
            sd = 1;
        end
        if ed > length(array_data(i).model_mean)
            ed = length(array_data(i).model_mean);
        end
        if p.fit_line
            t = array_data(i).model_T(sd:ed) - array_data(i).model_T(dex); 
            y = array_data(i).model_mean(sd:ed);
            output = polyfit(t',y,1);
            %output = fitlm(t',y,'Intercept', false);
            array_data(i).model_switch_to_0_strength(j) = output(1);
        else
            temp = mean(array_data(i).model_mean(dex:ed)) - mean(array_data(i).model_mean(sd:dex));
            array_data(i).model_switch_to_0_strength(j) = temp;
        end
    end

    for j=1:length(array_data(i).model_switch_to_1)
        dex = find(array_data(i).model_T == array_data(i).model_switch_to_1(j));   
        sd = dex - nw;
        ed = dex + nw;
        if sd < 1;
            sd = 1;
        end
        if ed > length(array_data(i).model_mean)
            ed = length(array_data(i).model_mean);
        end
        if p.fit_line
            t = array_data(i).model_T(sd:ed) - array_data(i).model_T(dex); 
            y = array_data(i).model_mean(sd:ed);
            output = polyfit(t',y,1);
            array_data(i).model_switch_to_1_strength(j) = output(1);
        else
            temp = mean(array_data(i).model_mean(dex:ed)) - mean(array_data(i).model_mean(sd:dex));
            array_data(i).model_switch_to_1_strength(j) = temp;
        end
    end
end
