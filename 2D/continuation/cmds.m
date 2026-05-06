%% demo for pBC 1D, clear workspace 
close all; keep pphome; 
%% cell 1: init
p=[]; par=[0.5 0.6 0]; % h_a, initial concentration, chemical potential (lagrange mult mass)
p=chinit(p,12.5,500,par); p.sw.qjac=0; %p.sw.verb=2;
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
% Solution vector
h = p.u(1:p.nu);
Q = h.^3 / 3;

x = getpte(p);
x = x(:);

fem = p.pdeo.fem;
gr  = p.pdeo.grid;

% --- compute logh
N = length(x);

slope = (h(2:end) - h(1:end-1)) ./ (x(2:end) - x(1:end-1));

dh_dx = zeros(N,1);
dh_dx(2:N-1) = 0.5*(slope(1:end-1) + slope(2:end));
dh_dx(1) = slope(1);
dh_dx(end) = slope(end);


% --- FEM assembly
[Kdiff,~,~] = fem.assema(gr, x.^2, 0, 0);
Kconv = convection(fem, gr, -x-3*x.^2.*dh_dx./h);
[~,Mreact,~] = fem.assema(gr, 0, 1, 0);

Krhs = convection(fem, gr, x.^2);

RHS = -Krhs * h;

LHS = -Kdiff + Kconv - Mreact;

psi = LHS \ RHS;

% Print residual
res = LHS*psi-RHS;
norm(res)
% Compute residual without FEM
dx = x(2)-x(1);

psix  = gradient(psi, dx);
psixx = gradient(psix, dx);

hx = gradient(h, dx);
dlogQdx = gradient(log(Q), dx);

R = psixx ...
    - dlogQdx .* psix ...
    + (1./x).*psix ...
    - (1./x.^2).*psi ...
    + hx;
norm(R,Inf)

figure(9);
plot(x,psi);
%p=cont(p,10); 
%branch = p.branch;
%writematrix(branch,'1D_c0_continuation_22_04_26.txt');