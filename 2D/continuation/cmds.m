%% demo for pBC 1D, clear workspace 
close all; keep pphome; 
%% cell 1: init
p=[]; par=[0.5 0.6 0]; % h_a, initial concentration, chemical potential (lagrange mult mass)
p=chinit(p,6.25,200,par); p.sw.qjac=0; %p.sw.verb=2;
%% Continuation parameters
p.nc.ilam = [2, 3];
p.nc.nq=1;
p.nc.lammax=2; p.sol.ds=0.005; p.nc.dsmax=0.01; 
%% First branch continuation
p=setfn(p,'tr'); p=findbif(p);
%% Switch to droplet branch
p.sol.ds=0.1;
p=swibra('tr','bpt1','b1',-0.1); p=cont(p,30); 
p.sw.bifcheck=0;
p.sw.foldcheck=1;

%%
% Solve ode for Psi v=1,2 using PDEtoolbox
h = p.u(1:p.nu); % Exclude the parameters of the pde
Qin = h.^3/3;
x=getpte(p); x=x';
% Derivatives
dhdr = p.mat.Kx*h; % or M/Kx*h ??
dQdr = dhdr.*(h.^2);
% I write the residual for the psi1 ode
% a(r)psi''+b(r)psi'+c(r)psi+d(r)=0
% As there is no non-linearity, I can solve psi by solving a linear system!
% (I can use just the \ operatioN!)
a_coeff = x.^2; b_coeff = x.*(1-x.*(dQdr./Qin)); c_coeff = -1*ones(p.nu,1); d_coeff = (x.^2).*dhdr;

% why multiply on right by x dependent coefficients?
LHS = -spdiags(a_coeff,0,p.np,p.np)*p.mat.K + spdiags(b_coeff,0,p.np,p.np)*p.mat.Kx - p.mat.M;
RHS = -d_coeff;
u = LHS \ RHS;
% Solve ode for psi1
%[rr, psi1] = ode45(@(r,psi) psi1_ode(r,psi,Qfun, hfun, Qrfun, hrfun), [p.vol,0], [0; 0]);

norm(LHS*u-RHS)


figure(9);
plot(x,u);
%p=cont(p,10); 
%branch = p.branch;
%writematrix(branch,'1D_c0_continuation_22_04_26.txt');