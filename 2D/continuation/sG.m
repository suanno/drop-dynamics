function r=sG(p,u)  % Alen Cahn with Natural bc plus lagrange multiplier for mass (chemical potential mu0)
par=u(p.nu+1:end); h=u(1:p.nu); % split u into parameters and PDE variables 
h_a=par(1); mu0=par(3);
x=getpte(p); x=x';
f=-x.*deriv_wetting_potential(h,h_a)+mu0*x;% Minus derivative of wetting potential + chemical potential (Lagrange multiplier of the mass)

F=p.mat.M*f;
% We consider the polar laplacian with radial symmetry (no theta
% derivative)
% The CM translation invariance is not necessary for radial equation as r=0 is
% a special point
r=p.mat.K*h-F;     % The convective derivative Kx has a + sign while the laplacian K has a - sign
  