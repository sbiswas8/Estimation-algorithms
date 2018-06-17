function PF_out = PFupdate(particles,wt,msmt,mdynamics,Ns)
n = length(particles(:,1));
PF_out.state = zeros(n,1);
PF_out.Pxx = zeros(n,n);
for i = 1:Ns
    X = particles(:,i);
    PF_out.state = PF_out.state + wt(i)*X;
end

for i = 1:Ns
    X = particles(:,i);
    PF_out.Pxx = PF_out.Pxx + wt(i)*(PF_out.state-X)*(PF_out.state-X)';
end
PF_out.Wt = wt;
PF_out.rsd = msmt.Z - mdynamics(PF_out.state);