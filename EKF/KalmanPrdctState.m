function state_p = KalmanPrdctState(state,dt,dynamics)
% Computes X^-
% called by of EKF function
% Author: Sanat K Biswas
% Australian Centre for Space Engineering Research
% UNSW Sydney
% December 2015: Last revision: 16-09-2017
% email: sanat@iiitd.ac.in
global t_step
tspan = [0 dt];
[T, S] = runge_kutta4(dynamics,tspan,state,t_step);
state_p = S(end,:)';