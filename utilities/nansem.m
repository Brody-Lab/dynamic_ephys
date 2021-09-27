function y=nansem(x)
%    nonanx=x(~isnan(x));
    y = nanstd(x) ./ sqrt(sum(~isnan(x)));
    
end
