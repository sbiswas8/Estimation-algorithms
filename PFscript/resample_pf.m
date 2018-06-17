function [prtcls,wt] = resample_pf(particles,weights,Ns)
% if sum(weights)~= 0
%     idx = randsample(1:Ns, Ns, true, weights);
% else
    edges = min([0 cumsum(weights)],1); % protect against accumulated round-off
    edges(end) = 1;                 % get the upper edge exact
    u1 = rand/Ns;
      % this works like the inverse of the empirical distribution and returns
      % the interval where the sample is to be found
    [~, idx] = histc(u1:1/Ns:1, edges);
%end
prtcls = particles(:,idx);
wt = repmat(1/Ns, 1, Ns);