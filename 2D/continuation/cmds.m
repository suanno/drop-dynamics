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
% Shift h(x) such that x=0 is the max
x = flip(x);
% Build matrices fo FEM
fem=p.pdeo.fem; gr=p.pdeo.grid; 
[K,~,~]=fem.assema(gr,1,1,1);
Kxr=convection(fem,gr,1/x);
[KdlogQdr,~,~]=fem.assema(gr,log(Qin),1,1);
[~,Mr2,~]=fem.assema(gr,1,1,1);
Kx=convection(fem,gr,1);
LHS = -K+Kxr+KdlogQdr-Mr2;
RHS = -Kx*h;
psi = LHS \ RHS;

% Print residual
res = LHS*psi-RHS;
norm(res)
% Compute residual without FEM
psix  = gradient(psi, x);        % first derivative
psixx = gradient(psix, x);       % second derivative
res = psixx + (1./x).*psix.*(1-x.*gradient(log(Qin),x)) - (1./x.^2).*psi + gradient(h,x);
norm(res)

figure(9);
plot(x,psi);
%p=cont(p,10); 
%branch = p.branch;
%writematrix(branch,'1D_c0_continuation_22_04_26.txt');