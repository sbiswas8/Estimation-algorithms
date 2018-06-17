function [particles,m] = CrtParticles(prtcls,dynamics,Ns,dt)
global t_step
%snoise = system_noise();
n = length(prtcls(:,1));
particles = zeros(n,Ns);
tspan = [0 dt];
m = 0;
for i = 1:Ns
    X = prtcls(:,i);
    [tout, Xpro] = runge_kutta4(dynamics,tspan,X,t_step);
    particles(:,i) = Xpro(end,:)';
    prcsd = i/Ns*100;
    if rem(i,10)==0
        msg= sprintf('%.0f percent particles created\n',prcsd);
        fprintf(repmat('\b',1,m));
        fprintf(msg);
        m=numel(msg);
    end
end