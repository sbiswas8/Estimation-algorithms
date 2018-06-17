function [H, rsd] = PrdctMsmtCrtH(state_p,msmt,msmtf,Mjacobian)
% Computes deltaZ (residue) and H matrix
% called by of EKF function
% Author: Sanat K Biswas
% Australian Centre for Space Engineering Research
% UNSW Sydney
% December 2015: Last revision: 16-09-2017
% email: sanat@iiitd.ac.in
Zp = msmtf(state_p);
H = Mjacobian(state_p);
rsd = msmt.Z - Zp;