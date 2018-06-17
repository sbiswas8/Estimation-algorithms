function [Pxz, Pvv, rsd] = UKFPrdctMsmt(X_priori,Xai_priori,msmt,msmtf)
% Computes deltaZ (residue), cross covariance matrix and measurement
% covariance matrix
% Author: Sanat K Biswas
% Australian Centre for Space Engineering Research
% UNSW Sydney
% December 2015: Last revision: 16-09-2017
% email: sanat@iiitd.ac.in

global n Wm0 Wm R
Wc0 =Wm0;
Wc = Wm;
%% Measurement prediction
Chai = nan(length(msmt.Z),2*n+1);
    
for count = 1:2*n+1
    Chai(:,count) = msmtf(Xai_priori(:,count));
end
Z = Wm0*Chai(:,1);
for count = 2:2*n+1
    Z = Z + Wm*Chai(:,count);
end

%% Cross covariance
Pzz = Wc0*(Chai(:,1) - Z)*(Chai(:,1) - Z)';
for count = 2:2*n+1
    Pzz = Pzz + Wc*(Chai(:,count) - Z)*(Chai(:,count) - Z)';
end
rsd = msmt.Z - Z;
Pvv = Pzz + R;

Pxz = Wc0*(Xai_priori(:,1) - X_priori)*(Chai(:,1) - Z)';
for count = 2:2*n+1
    Pxz = Pxz + Wc*(Xai_priori(:,count) - X_priori)*(Chai(:,count) - Z)';
end