function [particles,m] = ESPCrtParticles(state,Pxx,dynamics,FJacobian,Ns,dt)

global t_step
n = length(state);
particles = zeros(n,Ns);
tspan = [0 dt];
m = 0;
%% mean computation
X_mean = state;
%SigMat = repmat(X_mean,1,Ns) - prtcls;
SigMat = zeros(n,Ns);
for j = 1:n
    SigMat(j,:) = normrnd(0,sqrt(Pxx(j,j)),1,Ns);
end
%% mean propagation
[tout, Xpro] = runge_kutta4(dynamics,tspan,X_mean,t_step);
X_meanpro = Xpro(end,:)';
J = FJacobian(X_mean);
PHI = expm(J*dt);

for count = 1:Ns
    XaiN1 = X_meanpro + PHI*SigMat(:,count);
    delf = FJacobian(X_meanpro + SigMat(:,count)/2);
    delF = expm(delf);
    XaiN2 = X_meanpro + PHI*SigMat(:,count)/2 +...
            delF*SigMat(:,count)/2;
    particles(:,count) = 2*XaiN2 - XaiN1;
    
    prcsd = count/Ns*100;
    if rem(count,10)==0
        msg= sprintf('%.0f percent particles created\n',prcsd);
        fprintf(repmat('\b',1,m));
        fprintf(msg);
        m=numel(msg);
    end
end