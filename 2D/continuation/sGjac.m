function Gu=sGjac(p,u)  % PDE Jacobian for Allen Cahn with Natural bc
par=u(p.nu+1:end); u=u(1:p.nu); % params, and u on periodic domain 
h_a=par(1);
fu = -second_deriv_wetting_potential(u, h_a);
%fu=par(1)-3*u.^2; %Derivative of f respect to u in sG.m
Fu=p.mat.M*spdiags(fu,0,p.nu,p.nu);  % put derivatives into (sparse) matrix 

% Polar coordinate laplacian with radial symmetry (no theta derivative)
Gu=p.mat.K-p.mat.Kx-Fu;