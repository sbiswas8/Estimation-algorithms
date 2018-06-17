function [state_C, rsd] = EKF(state,msmt,model)
% EKF- Extended Kalman Filter algorithm
% This function takes a posteriori state vector and error covariance, next
% measurement vector
% and provides state estimate at the next time epoch.
% Syntax: [state_C, rsd] = EKF(state,msmt,model)
% Inputs:
%   state.X         -> a posteriori state vector
%   state.P         -> a posteriori error covariance
%   msmt.Z          -> measurement vector at the next time epoch
%   msmt.dt         -> measurement time interval (example: 1 sec, 2 sec etc)
%   model.dynamics  -> dynamics of the system (a vector differential 
%                      equation implementation i.e. X_dot =f(X,t))
%   model.Fjacobian -> Jacobian function of the system dynamics
%   model.msmtf     -> a vector measurment function Z = h(X)
%   model.Mjacobian -> Jacobian function (computes H matrix) for the
%                      measurement
% Outputs:
%   state_C.X -> state estimate
%   state_C.P -> error covariance
%   rsd       -> residue (= Z - Z^-)
% Other m-files required: Yes
%   KalmanPrdctState.m description: computes X^-
%   KalmanPrdctCov.m description: computes P^-
%   PrdctMsmtCrtH.m description: Computes deltaZ (residue) and H matrix
%   KalmanCrct.m description: Kalman Filter correction steps
% Note: Column vectors must be used as state and measurement
% Author: Sanat K Biswas
% Australian Centre for Space Engineering Research
% UNSW Sydney
% December 2015: Last revision: 16-09-2017
% email: sanat@iiitd.ac.in
dt = msmt.dt;
state_p = KalmanPrdctState(state.X,dt,model.dynamics);
P_p = KalmanPrdctCov(state.X,state.P,dt,model.Fjacobian);
[H, rsd] = PrdctMsmtCrtH(state_p,msmt,model.msmtf,model.Mjacobian);
state_C = KalmanCrct(state_p,P_p,rsd,H);