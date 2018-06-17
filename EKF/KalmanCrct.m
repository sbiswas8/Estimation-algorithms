function [state_c] = KalmanCrct(state_p,P_p,rsd,H)
% Kalman Filter correction steps
% called by of EKF function
% Author: Sanat K Biswas
% Australian Centre for Space Engineering Research
% UNSW Sydney
% December 2015: Last revision: 16-09-2017
% email: sanat@iiitd.ac.in
global R
K = P_p*H'*(R + H*P_p*H')^-1;
state_c.X = state_p + K*rsd;
state_c.P = (eye(length(state_c.X)) - K*H)*P_p*(eye(length(state_c.X)) - K*H)' + K*R*K';
state_c.K = K;
for i = 1:length(state_c.X)
    state_c.STD(i) = sqrt(state_c.P(i,i));
end