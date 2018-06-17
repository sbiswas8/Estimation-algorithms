function state_c = SPUKFCrct(X_priori,P_p,rsd,Pxz,Pvv)
% SPUKF correction steps
% Author: Sanat K Biswas
% Australian Centre for Space Engineering Research
% UNSW Sydney
% December 2015: Last revision: 16-09-2017
% email: sanat@iiitd.ac.in
%% Gain calculation and update
K = Pxz*Pvv^-1;
X_post = X_priori + K*rsd;
P = P_p - K*Pvv*K';
% P = 0.5*P + 0.5*P';
% P = P + 1e-8*eye(length(P));
state_c.X = X_post;
state_c.P = P;
for i = 1:length(state_c.X)
    state_c.STD(i) = sqrt(state_c.P(i,i));
end
