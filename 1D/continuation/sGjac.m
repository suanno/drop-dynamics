function Gu=sGjac(p,u)  % PDE Jacobian for AC with pBC 
Hp=u(1:p.nu); % params, and u on periodic domain 
H=p.mat.fill*Hp; % extend ('fill') u to full domain 

fu = second_deriv_wetting_potential(H);
%fu=par(1)-3*u.^2; %Derivative of f respect to u in sG.m
Fu=p.mat.M0*(spdiags(fu,0,p.np,p.np)*p.mat.fill); % map fu to per.dom
%Gu=p.mat.K+par(4)*p.mat.Kx-Fu; % Translational invariance
Gu=p.mat.K+Fu;  % Center of mass at x=0