function [X_priori, Xai_priori] = ESPUKFPrdctState(Xai,dt,dynamics,FJacobian)
% computes X^- using Extrapolated Single Propagation Technique
% Author: Sanat K Biswas
% Australian Centre for Space Engineering Research
% UNSW Sydney
% December 2015: Last revision: 16-09-2017
% email: sanat@iiitd.ac.in
global t_step m n Wm0 Wm
tspan = [0 dt];
%% Propagation of Sigma points
Xai_priori = nan(n/2,2*n+1);
[tout, Xpro] = runge_kutta4(dynamics,tspan,Xai.X,t_step);
Xai_priori(:,1) = Xpro(end,1:m)';
J = FJacobian(Xai.X(1:m,1));
PHI = expm(J*dt);
for count = 2:n+1
    XaiN1 = Xai_priori(:,1) + PHI*Xai.SigMat(1:m,count-1);
    delf = FJacobian(Xai.X(1:m,1) + Xai.SigMat(1:m,count -1)/2);
    delF = expm(delf);
    XaiN2 = Xai_priori(:,1) + PHI*Xai.SigMat(1:m,count-1)/2 +...
            delF*Xai.SigMat(1:m,count-1)/2;
    Xai_priori(:,count) = 2*XaiN2 - XaiN1 + Xai.SigMat(m+1:n,count - 1)*dt;
end

for count = n+2: 2*n+1
    XaiN1 = Xai_priori(:,1) - PHI*Xai.SigMat(1:m,count-1-n);
    delf = FJacobian(Xai.X(1:m,1) - Xai.SigMat(1:m,count -1-n)/2);
    delF = expm(delf);
    XaiN2 = Xai_priori(:,1) - PHI*Xai.SigMat(1:m,count-1-n)/2 -...
            delF*Xai.SigMat(1:m,count-1-n)/2;
    Xai_priori(:,count) = 2*XaiN2 - XaiN1 - Xai.SigMat(m+1:n,count - 1-n)*dt;
end
    %UKFsigpt = Xai_priori;
X_priori = Wm0*Xai_priori(:,1);
for count = 2:2*n+1
    X_priori = X_priori + Wm*Xai_priori(:,count);
end