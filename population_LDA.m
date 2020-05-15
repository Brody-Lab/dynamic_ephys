function population_LDA(which_switch, savefig) 
%%%%%% LOOK AT ME, READ ME FIRST. This approach isn't going to work because the method requires knowing the covariance between each cell. we don't know that, and can't know that for cells that are not recorded simulatneously. not sure how to proceed. 

%which_switch = 'generative';
if nargin < 2
    savefig = 0;
end
[res, cellids, computed,not_computed] = load_LDA_population(which_switch, 0);
lags = res{1}.lags > -.5 & res{1}.lags < 1;
% iterate over time points
%%for i = 1:sum(lags)
i = 100;
    %at each time point, make x and T vectors
    fr1 = [];
    fr2 = [];
    for j=1:sum(computed)
         fr1= [fr1; res{j}.STR_right_real(:,i)];
         fr2= [fr2; res{j}.STR_left_real(:,i)];
    end
    fr1(isnan(fr1)) = [];
    fr2(isnan(fr2)) = [];
    fr = [fr1; fr2];
    states = [zeros(size(fr1)); ones(size(fr2))];
    % do LDA

    [Y,W,LAMBDA] = LDA(fr,states);
    % compute Fisher score
    % save


end



