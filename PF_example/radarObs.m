function z = radarObs(X)
M = 1e5;
H = 1e5;
z = sqrt(M^2 +(X(1) - H).^2);