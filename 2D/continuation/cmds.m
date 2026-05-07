%% demo for pBC 1D, clear workspace 
close all; keep pphome; 
%% cell 1: init
p=[]; par=[0.5 0 0]; % h_a, initial concentration, chemical potential (lagrange mult mass)
p=chinit(p,12.5,2000,par); p.sw.qjac=0; %p.sw.verb=2;
%% Continuation parameters
p.nc.ilam = [2, 3];
p.nc.nq=1;
p.nc.lammax=12.5; p.sol.ds=0.01; p.nc.dsmax=0.5;
%% First branch continuation
p=setfn(p,'tr'); p=findbif(p,30);
%% Switch to droplet branch
p.sol.ds=+0.01;
p=swibra('tr','bpt1','b1',+0.01); p=cont(p,30); 
p.sw.bifcheck=0;
p.sw.foldcheck=1;

%%


% Solve ode for Psi v=1,2 using PDEtoolbox
h = p.u(1:p.nu); % Exclude the parameters of the pde
Qin = h.^3/3;
x=getpte(p); x=x';

% Compute residuals for h0in ODE computing h0out from hmax
h_a=par(1);
% hr=gradient(h,x);
% hrr=gradient(hr,x);
dwout=(deriv_wetting_potential(max(h),h_a)-deriv_wetting_potential(h_a,h_a))/(max(h)-h_a);
% res=-x.*hrr-hr+x.*deriv_wetting_potential(h,h_a)-x.*dwout;
% norm(res)
fem=p.pdeo.fem; gr=p.pdeo.grid; 
[Kr,~,~]=fem.assema(gr,x,1,1);
Kx=convection(fem,gr,1);
[~,Mr,~]=fem.assema(gr,1,1,1);
res=(Kr-Kx)*h+Mr*(deriv_wetting_potential(h,h_a)-dwout);
norm(res)


drlogQin=3*(gradient(h,x)./h);
% Build matrices for FEM
fem=p.pdeo.fem; gr=p.pdeo.grid; 
[Kr2,~,~]=fem.assema(gr,x.^2,1,1);
Kxr2=convection(fem,gr,x.^2);
Kx2r=convection(fem,gr,2*x);
Kxr2logQin=convection(fem,gr,x.*(1-x.*drlogQin));
[~,M,~]=fem.assema(gr,1,1,1);
LHS = Kxr2logQin--Kx2r-Kr2-M;
RHS = -Kxr2*h;
psi = LHS \ RHS;

% Print residual
res = LHS*psi-RHS;
norm(res)
% Compute residual without FEM
% psix  = gradient(psi, x);        % first derivative
% psixx = gradient(psix, x);       % second derivative
% res = (x.^2).*psixx + x.*psix.*(1-x.*drlogQin) - psi + (x.^2).*gradient(h,x);
% norm(res)

figure(9);
plot(x,psi);
%p=cont(p,10); 
%branch = p.branch;
%writematrix(branch,'1D_c0_continuation_22_04_26.txt');