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
p=swibra('tr','bpt1','b1',-0.1); p=cont(p,3); 
p.sw.bifcheck=0;
p.sw.foldcheck=1;
p=cont(p,500); 
branch = p.branch;
writematrix(branch,'1D_c0_continuation_22_04_26.txt');