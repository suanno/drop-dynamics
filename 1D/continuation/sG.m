function r=sG(p,u)  % AC with periodic BC plus lagrange multiplier for mass (chemical potential mu0)
par=u(p.nu+1:end); up=u(1:p.nu); % params, and u on periodic domain 
u=p.mat.fill*up; % extend ('fill') u to full domain 
h_a=par(1); mu0=par(3);
f=-deriv_wetting_potential(u,h_a)+mu0;% Minus derivative of wetting potential + chemical potential (Lagrange multiplier of the mass)
%f=lambda*u-u.^3+mu0;  par(2)>0
F=p.mat.M0*f; % multiply by M, map back to active nodes of periodic domain 
%r=p.mat.K*up+par(4)*p.mat.Kx*up-F; % Translation invariance
x=getpte(p); x=x'; % extract the point coordinates from p
r=p.mat.K*up-F-p.mat.M0*x*par(4);     % Fix center of mass at x=0
  