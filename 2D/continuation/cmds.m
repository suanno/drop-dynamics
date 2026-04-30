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

% Solve ode for Psi v=1,2 using PDEtoolbox
h = p.u(1:p.nu); % Exclude the parameters of the pde
Qin = h.^3;
x=getpte(p); x=x';
% Interpolate
dhdr = p.mat.M \ (spdiags(1/x,0,p.np,p.np)*p.mat.Kx*h);
hfun = @(r) interp1(x,h,r,'pchip',0);
hrfun = @(r) interp1(x,dhdr,r,'pchip',0);
Qfun = @(r) interp1(x,Qin,r,'pchip',0);
Qrfun = @(r) hrfun(r)*hfun(r)^2/3;
% Solve ode for psi1
[rr, psi1] = ode45(@(r,psi) psi1_ode(r,psi,Qfun, hfun, Qrfun, hrfun), [p.vol,0], [0; 0]);

figure(9);
plot(rr,psi1(:,1));
%p=cont(p,10); 
%branch = p.branch;
%writematrix(branch,'1D_c0_continuation_22_04_26.txt');