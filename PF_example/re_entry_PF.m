function [PFdata,pf_msgl] = re_entry_PF(Filter_parameter)
dbstop if 'error'
%load Filter_parameter
load Re_entry_ref
%% Filter Initialization
X_ini = [3e5; 2e4; 3e-5];
%P_ini = [1e6 0 0; 0 4e6 0; 0 0 1e-4];
%P_ini = [0.1 0 0; 0 0.01 0; 0 0 1e-4];
P_ini = Filter_parameter.PF.P_ini;
n = length(X_ini);
X_est(1,:) = X_ini';
cov(1,:) = sqrt([P_ini(1,1) P_ini(3,3) P_ini(3,3)]);
%% Model declaration
model.sys = @reentryDyn;
model.sysjacobian = @jacobian;
model.meas = @radarObs;
model.Pzx = @re_entryP_zx;
%model.sysnoise = @re_entry_sysNoise;

%% Global variables
global t_step Qk R
t_step = 0.1;
% w1 = 1e-8;
% w2 = 1e-8;
% w3 = 1e-8;
Qk = Filter_parameter.PF.Q;
R = Filter_parameter.PF.R;

%% estimation
Ns = Filter_parameter.Ns;
Sample_Info.Ns = Ns;
Sample_Info.thld = 0.5;
RSD = [];
m = 0;
k = 1;
prtcls = struct([]);
profile on
for i = 1:length(TIME) - 1
    if i == 1
        PF_in.state = X_ini;
        PF_in.Pxx = P_ini;
        PF_in.Wt = repmat(1/Ns,1,Ns);
        Xai = zeros(n,Ns);
        for j = 1:n
            Xai(j,:) = normrnd(0,sqrt(P_ini(j,j)),1,Ns);
        end
        PF_in.particles = repmat(X_ini,1,Ns) + Xai;
    else
        PF_in.state = PF_out.state;
        PF_in.Pxx = PF_out.Pxx;
        PF_in.Wt = PF_out.Wt;
        PF_in.particles = PF_out.particles;
    end
    msmt.Z = Z(i+1);
    msmt.dt = dt;
    [PF_out,msgl,N_eff] = particle_filter(PF_in,msmt,model,Sample_Info);
    X_est(i+1,:) = PF_out.state';
    if rem(i,10) == 0
        prtcls(k).pt = PF_out.particles;
        prtcls(k).ptnr = PF_out.prtcls_n_rsmpl;
        k = k+1;
    end
    cov(i+1,:) = sqrt([PF_out.Pxx(1,1) PF_out.Pxx(2,2) PF_out.Pxx(3,3)]);
    RSD = [RSD PF_out.rsd];
    
    prcsd = i/length(TIME)*100;
    msg= sprintf('Particle Filter: %.00f percent processed\n residue: %f\n sampling quality: %f\n ',prcsd,PF_out.rsd,N_eff);
    fprintf(repmat('\b',1,m+msgl));
    fprintf(msg);
    m=numel(msg);
end
pf_msgl = m;
PR = profile('info');
PFdata.ex_time = PR.FunctionTable(6).TotalTime/PR.FunctionTable(6).NumCalls;
PFdata.error = X_true - X_est;
PFdata.prtcls = prtcls;
PFdata.rsd = RSD;
%save(file_name,'PFdata')