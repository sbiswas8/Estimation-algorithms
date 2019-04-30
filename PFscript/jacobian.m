function delF = jacobian(Xai)
lambda = 5e-5;
delF = [0 -1 0; lambda*exp(-lambda*Xai(1,1))*Xai(2,1)^2*Xai(3,1)...
        -2*exp(-lambda*Xai(1,1))*Xai(2,1)*Xai(3,1) -exp(-lambda*Xai(1,1))*Xai(2,1)^2; 0 0 0];