function [hit] = compute_linear_agent(trial,p,lambda)

lam = lambda; 
dt = p.optimal.dt;
y = 0;
y = zeros(round(trial.T/dt),1);
L = round(trial.left./dt);
R = round(trial.right./dt);
ev = zeros(round(trial.T/dt),1);

L(L==0) = 1;
R(R==0) = 1;

L(L>length(ev)) = length(ev);
R(R>length(ev)) = length(ev);

ev(R) = 1;
ev(L) = ev(L) -1;

% Iterate trial
for t=1:round(trial.T/dt)-1;
    y(t+1) = y(t) + ev(t) + lam*dt*y(t);
end
hit = y(end)>0;

