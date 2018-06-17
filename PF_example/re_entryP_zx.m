function p = re_entryP_zx(msmt,X,Pxx,mdynamics)
global R

M = 1e5;
H = 1e5;

z_cap = mdynamics(X);
Hx = [-2*X(1)/(sqrt(M^2 +(X(1) - H).^2)) 0 0];
Pzx = Hx*Pxx*Hx'+R;
p = normpdf(msmt.Z-z_cap,0,sqrt(R));