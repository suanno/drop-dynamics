function Gu=sGjac(p,u)  % PDE Jacobian for Allen Cahn with Natural bc
H=u(1:p.nu); % params, and u on periodic domain 
x=getpte(p); x=x';
fu = -x.*second_deriv_wetting_potential(H);
%fu=par(1)-3*u.^2; %Derivative of f respect to u in sG.m
Fu=p.mat.M*spdiags(fu,0,p.nu,p.nu);  % put derivatives into (sparse) matrix 


% Polar coordinate laplacian with radial symmetry (no theta derivative)
Gu=p.mat.K-Fu;