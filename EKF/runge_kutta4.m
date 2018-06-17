function [T, X] = runge_kutta4(func,t_span,R,t_step)
%Simple RK-4 implementation
% Author: Sanat K Biswas
% Australian Centre for Space Engineering Research
% UNSW Sydney
% December 2015: Last revision: 16-09-2017
% email: sanat@iiitd.ac.in
h = t_step;
t_sim = t_span(1):h:t_span(2);
X = nan(length(t_sim),length(R));
X(1,:) = R';
for i=1:(length(t_sim))-1   %calculation loop
    k_1 = func(t_sim(i),R);
    k_2 = func(t_sim(i)+0.5*h,R +0.5*h*k_1);
    k_3 = func((t_sim(i)+0.5*h),(R+0.5*h*k_2));
    k_4 = func((t_sim(i)+ h),(R +k_3*h));

    R = R + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;  %main equation
    X(i+1,:) = R';
end
T = t_sim;