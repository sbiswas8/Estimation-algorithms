function Xai = SPUKFCrtDelSig(X,Px)
% Computes delta sigma points
% Author: Sanat K Biswas
% Australian Centre for Space Engineering Research
% UNSW Sydney
% December 2015: Last revision: 16-09-2017
% email: sanat@iiitd.ac.in
global n lambda_kf
%% Delta Sigma point calculation
Xai.SigMat  = sqrtm((n + lambda_kf)*Px);
Xai.X = X;