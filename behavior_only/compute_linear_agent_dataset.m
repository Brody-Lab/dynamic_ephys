function [S] = compute_linear_agent_dataset(rat,S,p);

if ~isfield(S.pd{1}.bupsdata{1},'noiseAnswer')
lambda = p.optimal.lambda;
for i=1:length(S.pd);
    for j=1:length(S.pd{i}.hits)
        S.pd{i}.bupsdata{j}.noiseAnswer = compute_linear_agent(S.pd{i}.bupsdata{j},p,lambda);
    end
end

%save(['rat_data/Haz_1/model_data/' rat '_full_dataset.mat'],'S','p','rat');
end

