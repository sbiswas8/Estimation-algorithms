function [weights, N_eff] = CalcWeights(wt,Pxx,particles,msmt,mdynamics,Pzxmodel,Ns)
weights = zeros(1,Ns);
for i = 1:Ns
    X = particles(:,i);
    Pzx = Pzxmodel(msmt,X,Pxx,mdynamics);
    weights(i) = wt(i)*Pzx;
end
Swt = sum(weights);
%% normalization
for i = 1:Ns
    if Swt~=0
        weights(i) = weights(i)/Swt;
    end
end
N_eff = 1/sum(weights.^2);