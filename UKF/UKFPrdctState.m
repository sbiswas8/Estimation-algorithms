function [X_priori, Xai_priori] = UKFPrdctState(Xai,dt,dynamics)
% computes X^-
% Author: Sanat K Biswas
% Australian Centre for Space Engineering Research
% UNSW Sydney
% December 2015: Last revision: 16-09-2017
% email: sanat@iiitd.ac.in

global t_step m n Wm0 Wm
tspan = [0 dt];
%% Propagation of Sigma points
Xai_priori = nan(n/2,2*n+1);
for count = 1: 2*n+1
    [tout, Xpro] = runge_kutta4(dynamics,tspan,Xai(:,count),t_step);
    Xai_priori(:,count) = Xpro(end,1:m)';%rem(Xpro(end,1:m)',1*pi);
end
    %UKFsigpt = Xai_priori;
X_priori = Wm0*Xai_priori(:,1);
for count = 2:2*n+1
    X_priori = X_priori + Wm*Xai_priori(:,count);
end