function p=chinit(p,lx,nx,par)  % init routine for AC on interval with pBC 
p=stanparam(p); screenlayout(p); p.sw.sfem=-1; 
p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; p.fuha.qf=@qf;
%p.fuha.e2rs=@e2rs; % the relevant fun.handles
pde=stanpdeo1D(lx,2*lx/nx); p.pdeo=pde; % domain and mesh
p.np=pde.grid.nPoints; p.nu=p.np; p.sol.xi=1/(p.nu); 
%[po,t,e]=getpte(p);
%p.mesh.bp=po; p.mesh.bt=t; p.mesh.be=e; % background mesh (for mesh adaption) 
c0=par(2);
p.u=c0*ones(p.np,1); p.u=[p.u; par']; % initial guess (homogeneous) 
p.vol=2*lx;   % Total volume
p.sw.bcper=1; p=box2per(p); % prepare fill, drop for periodic BC, here in x
p.nc.nsteps=20; p.sw.foldcheck=1; p.plot.auxdict={'lambda','c0','mu0'}; 
p.plot.pstyle=1; p.usrlam=[0 0.5 1]; p.nc.nsteps=100; p.sw.jac=1; 
p.sw.bifcheck=2;
p.fuha.outfu = @bra;      % Measured observables along the branch
p.spcontsw=0;   % Compute numerically spjac.m 