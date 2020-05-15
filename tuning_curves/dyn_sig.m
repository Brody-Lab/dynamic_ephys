% y=sig4(beta,x)
% y0=beta(1);
% a=beta(2);
% x0=beta(3);
% b=beta(4);

function y=dyn_sig(beta,x)

%y0=beta(1);
%a=beta(2);
%x0=beta(3);
%b=1./beta(4);
x0 = beta(1);
b = 1./beta(2);

y=1./(1+ exp(-(x-x0)./b));
% slope at x = 0 given by beta(2)/4
