function P_p = ESPUKFPrdctCov(X_priori,Xai_priori)
% Computes P^-
% Author: Sanat K Biswas
% Australian Centre for Space Engineering Research
% UNSW Sydney
% December 2015: Last revision: 16-09-2017
% email: sanat@iiitd.ac.in
global n Wm0 Wm
Wc0 = Wm0;% + (1 - alpha^2 + beta);
Wc = Wm;

P_p = Wc0*(Xai_priori(:,1) - X_priori)*(Xai_priori(:,1) - X_priori)';

for count = 2:2*n+1
    P_p = P_p + Wc*(Xai_priori(:,count) - X_priori)*(Xai_priori(:,count) - X_priori)';
end