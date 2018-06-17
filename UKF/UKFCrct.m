function state_c = UKFCrct(X_priori,P_p,rsd,Pxz,Pvv)
% UKF correction steps
% Author: Sanat K Biswas
% Australian Centre for Space Engineering Research
% UNSW Sydney
% December 2015: Last revision: 16-09-2017
% email: sanat@iiitd.ac.in

%% Gain calculation and update
global det_K
K = Pxz*Pvv^-1;
X_post = X_priori + K*rsd;
P = P_p - K*Pvv*K';
state_c.X = X_post;%rem(X_post,19*pi);
state_c.P = P;
state_c.K =K;
for i = 1:length(state_c.X)
    state_c.STD(i) = sqrt(state_c.P(i,i));
end
det_K = K;
