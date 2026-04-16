%% demo for pBC 1D, clear workspace 
close all; keep pphome; 
%% cell 1: init
p=[]; par=[0.5 0.62 0 0]; % h_a, initial concentration, chemical potential (lagrange mult mass), velocity (lagrange mult translational invariance OR center of mass lagrange mult)
p.reducedmass=0;
p.xcm = 0;
p=chinit(p,25,200,par); p.nc.nq=1; p.sw.qjac=0; %p.sw.verb=2;
%% Continuation parameters
p.nc.ilam = [2 3];
p.nc.lammax=0.64; p.sol.ds=0.001; p.nc.dsmax=0.01; 
%% First branch continuation
p=setfn(p,'tr'); p=findbif(p);
%% Switch to 1 droplet branch
p.sol.ds=-0.01;
p.nc.dsmin=-0.1;
p.nc.dsmax=-0.001;
p=swibra('tr','bpt1','b1',-0.1); p=cont(p,1); 
p.sw.bifcheck=0;
p.sw.foldcheck=1;
%% Add the translational invariance AND go towards a small droplet compared to the system size
% Compute the CM position
x=getpte(p); x=x'; % extract the point coordinates from p
up = p.u(1:p.nu);
u=p.mat.fill*up;
p.xcm = sum(p.mat.M0*(u.*x(1:end)));
%p.sol.ds=0.01;
p.nc.nq=2; p.nc.ilam=[2 3 4]; p.tau=[p.tau; 0];
p=cont(p,30);
%%
p=swibra('b1','fpt1','f1',-0.1); %Negative upper branch positive lower branch
p.nc.lammax=1;
p.sw.bifcheck=0;
p.sw.foldcheck=0;
p.nc.lammax=5;
p=cont(p,500);
%% [NOT NECESSARY, I can measure the lagrange multiplier while continuating in mass!]
% Remove mass constraint (but not xcm) and continuation in lagrange
% multiplier of the mass
%p.nc.nq=1;
%p.fuha.qf=@qf_nomass;
%p.nc.ilam=[3 4];
%p.nc.lammax=5;
%p=cont(p,30);



