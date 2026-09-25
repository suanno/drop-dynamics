%% demo for pBC 1D, clear workspace 
close all; keep pphome; 
%% cell 1: init
p=[]; par=[1 0 0]; % concentration, chemical potential (lagrange mult mass), velocity (lagrange mult translational invariance OR center of mass lagrange mult)
p.reducedmass=0;
p.xcm = 0;
p=chinit(p,25,500,par); p.nc.nq=1; p.sw.qjac=0; %p.sw.verb=2;
%% Continuation parameters
p.nc.ilam = [1 2];
p.nc.lammax=2; p.sol.ds=0.005; p.nc.dsmax=0.01; 
%% First branch continuation
p=setfn(p,'tr'); p=findbif(p);
%% Switch to 1 droplet branch
%p.sol.ds=-0.01;
%p.nc.dsmin=-0.1;
%p.nc.dsmax=-0.001;
p=swibra('tr','bpt1','b1',-0.1); p=cont(p,1); 
p.sw.bifcheck=0;
p.sw.foldcheck=1;
p=cont(p,400);
writematrix(p.branch,'1D_c0_continuation_nbc.txt');
%% The conservation of Xcm is not necessary for nbc bc as one of the boundary is the max of the drop