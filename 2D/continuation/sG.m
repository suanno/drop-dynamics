function r=sG(p,u)  % Alen Cahn with Natural bc plus lagrange multiplier for mass (chemical potential mu0)
par=u(p.nu+1:end); u=u(1:p.nu); % split u into parameters and PDE variables 
h_a=par(1); mu0=par(3);
f=-deriv_wetting_potential(u,h_a)+mu0;% Minus derivative of wetting potential + chemical potential (Lagrange multiplier of the mass)
x=getpte(p); x=x';

F=p.mat.M*f;
% We consider the polar laplacian with radial symmetry (no theta
% derivative)
% The CM translation invariance is not necessary for radial equation as r=0 is
% a special point
r=p.mat.K*u-p.mat.Kx*u-F;     % The convective derivative Kx has a + sign while the laplacian K has a - sign
  