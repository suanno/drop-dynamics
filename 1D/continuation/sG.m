function r=sG(p,u)  % AC with periodic BC plus lagrange multiplier for mass (chemical potential mu0)
par=u(p.nu+1:end); Hp=u(1:p.nu); % params, and u on periodic domain 
H=p.mat.fill*Hp; % extend ('fill') u to full domain 

f=deriv_wetting_potential(H)-par(2);% Minus derivative of wetting potential + chemical potential (Lagrange multiplier of the mass)
F=p.mat.M*f; % multiply by M, map back to active nodes of periodic domain 

x=getpte(p); x=x'; % extract the point coordinates from p
r=p.mat.K*H+F+p.mat.M*x*par(3);     % Fix center of mass at x=0
  