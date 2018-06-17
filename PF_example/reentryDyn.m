function x_dot = reentryDyn(t,x)
global Qk
%% Parameters
lambda = 5e-5;

%% equations

x_dot(1,1) = -x(2) + normrnd(0,sqrt(Qk(1,1)));
x_dot(2,1) = -exp(-lambda*x(1))*x(2)^2*x(3) + normrnd(0,sqrt(Qk(2,2)));
x_dot(3,1) = normrnd(0,sqrt(Qk(3,3)));
%x_dot(3,1) = 0;