function P_p = KalmanPrdctCov(state,P,dt,Fjacobian)
% Computes P^-
% called by of EKF function
% Author: Sanat K Biswas
% Australian Centre for Space Engineering Research
% UNSW Sydney
% December 2015: Last revision: 16-09-2017
% email: sanat@iiitd.ac.in
global Qk
J = Fjacobian(state);
PHI = expm(J*dt);
P_p = PHI*P*PHI' + PHI*Qk*PHI'*dt;