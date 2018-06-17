function Xai = CrtSigmaPt(X,Px)
% Computes sigma points
% Author: Sanat K Biswas
% Australian Centre for Space Engineering Research
% UNSW Sydney
% December 2015: Last revision: 16-09-2017
% email: sanat@iiitd.ac.in
global n lambda_kf
%% Sigma point calculation
SigMat  = sqrtm((n + lambda_kf)*Px);
Xai = nan(n,2*n+1);
Xai(:,1) = X;
for count = 2:n+1
    Xai(:,count) = X + SigMat(:,count - 1);
end
for count = n+2:2*n+1
    Xai(:,count) = X - SigMat(:,count - 1 - n);
end