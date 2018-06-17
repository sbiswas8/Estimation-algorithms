function [state_C, rsd] = UKF(state,msmt,model)
% UKF- Unscented Kalman Filter algorithm
% This function takes a posteriori state vector and error covariance, next
% measurement vector and provides state estimate at the next time epoch 
% using UKF.
% Syntax: [state_C, rsd] = UKF(state,msmt,model)
% Inputs:
%   state.X        -> a posteriori state vector
%   state.P        -> a posteriori error covariance
%   msmt.Z         -> measurement vector at the next time epoch
%   msmt.dt        -> measurement time interval (example: 1 sec, 2 sec etc)
%   model.dynamics -> dynamics of the system (a vector differential 
%                     equation implementation i.e. X_dot =f(X,t))
%   model.msmtf    -> a vector measurment function Z = h(X)
% Global variables: Qk (process noise matrix) must be defined in the main script
%                   rest of the global variables are not accessible from
%                   the main script unless declared there
% Outputs:
%   state_C.X -> state estimate
%   state_C.P -> error covariance
%   rsd       -> residue (= Z-Z^-)
% Other m-files required: Yes
%   CrtSigmaPt.m description: computes sigma points
%   UKFPrdctState.m description: computes X^-
%   UKFPrdctCov.m description: computes P^-
%   UKFPrdctMsmt.m description: Computes deltaZ (residue)
%   UKFCrct.m description: UKF correction steps
% Note: Column vectors must be used as state and measurement
% Author: Sanat K Biswas
% Australian Centre for Space Engineering Research
% UNSW Sydney
% December 2015: Last revision: 16-09-2017
% email: sanat@iiitd.ac.in
global Qk m n Wm0 Wm lambda_kf
m = length(state.X);
n = 2*m;
%% Parameters
alpha = sqrt(3/n);
beta = 2;
Kai = 0;
lambda_kf = alpha^2*(n + Kai) - n;
Wm0 = lambda_kf/(n + lambda_kf);
Wm = 1/(2*(n  + lambda_kf));
%% Augmentation
X = [state.X;zeros(m,1)];
Px = [state.P zeros(m);...
     zeros(m) Qk];
%% Filter
dt = msmt.dt;
Xai = CrtSigmaPt(X,Px);
[X_priori, Xai_priori] = UKFPrdctState(Xai,dt,model.dynamics);
P_p = UKFPrdctCov(X_priori,Xai_priori);
[P_XZ, Pvv, rsd] = UKFPrdctMsmt(X_priori,Xai_priori,msmt,model.msmtf);
state_C = UKFCrct(X_priori,P_p,rsd,P_XZ,Pvv);