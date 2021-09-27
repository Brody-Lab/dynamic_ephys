function [data_auc, sd_p, ci]=bootroc(A,B,BOOTS, CI)
% [auc, auc_pctile, shuff_ci]=bootroc(A,B,BOOTS, CI);
%
% A is a vector of the value of elements from condition A
% B is a vector of the value of elements from condition B
% 

% Modified Aug 13 by tyler so that the pvalue is not the probability that
% the real auc came from the shuffle distribution, but the percentile of the
% auc value in the shuffle distribution (the old method is highly
% dependent on number of shuffles). Also, relabeling the third output as
% shuff_ci rather than auc ci, to clarify that it represents the confidence

% modified 1/20/20 by tyler so that names are more informative

% intervals around the shuffle distribution.
if nargin<3
	BOOTS=1000;
end

if nargin<4
  CI=95;
  alpha = 100-CI;
end


A=A(~isnan(A));
B=B(~isnan(B));

if isempty(A) || isempty(B)
    data_auc=nan;
    sd_p=nan;
    boot_auc=nan;
else

data_auc = auc(A,B);
sA=numel(A);
ALL_DATA=[A(:);B(:)];
boot_auc=0.5+zeros(BOOTS,1);
parfor bx=1:BOOTS

	shuff_d=ALL_DATA(randperm(numel(ALL_DATA)));
	A=shuff_d(1:sA);
	B=shuff_d(sA+1:end);

	boot_auc(bx)=auc(A,B);
end

%sd_p=get_p(sd, boot_score);
sd_p = mean(boot_auc <= data_auc);
end

ci=prctile(boot_auc, [alpha/2 (100-alpha/2)]);