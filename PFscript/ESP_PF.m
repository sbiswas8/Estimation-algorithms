function [PF_out,msgl,N_eff] = ESP_PF(PF_in,msmt,model,Sample_Info)
%% Description
%% Inputs: 
% PF_in structure
% PF_in.state = previous posteriori estimate
% PF_in.Pxx = previous posteriori covariance
% PF_in.Wt = previous particle weights
% PF_ini.particles = previous particles

% model structure
% model.sys = system model
% model.meas = measurement model
% model.sys_noise = system noise generation
% model.Pzx = function to create condtional probality P of Z given x

% msmt structure
% msmt.Z = observation
% msmt.dt = measurement interval

% Sample_Info structure
% Sample_Info.Ns = no. of samples
% Sample_Info.thld = threshold

%% Output:
% PF_out.state = estimated state
% PF_out.Pxx = estimated covariance
% PF_out.Wt = calculated weights
% PF_out.rsd = residue
%%
Ns = Sample_Info.Ns;
Nt = Sample_Info.thld*Ns;
dt = msmt.dt;
[particles,msgl] = ESPCrtParticles(PF_in.state,PF_in.Pxx,model.sys,model.sysjacobian,Ns,dt);
[weights, N_eff] = CalcWeights(PF_in.Wt,PF_in.Pxx,particles,msmt,model.meas,model.Pzx,Ns);
prtcls_no_rsmpl = particles;
if N_eff < Nt || N_eff == inf
    [particles,weights] = resample_pf(particles,weights,Ns);
end
PF_out = PFupdate(particles,weights,msmt,model.meas,Ns);
PF_out.prtcls_n_rsmpl = prtcls_no_rsmpl;
PF_out.particles = particles;