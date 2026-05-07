function p=chinit(p,lx,nx,par)  % init routine for AC on interval with pBC 
p=stanparam(p); screenlayout(p); p.sw.sfem=-1; 
p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; p.fuha.qf=@qf;
% We need the grid to be [0, 2*Lx] not [-Lx, Lx] as the x-coordinate is the
% radial coordinate for polar coordinates
pde=stanpdeo1D(lx,2*lx/nx); p.pdeo=pde; % symmetric [-Lx,Lx] domain and mesh
pde.grid.interval([0, 2*lx], 2*lx/nx); % asymmetric [0, 2*Lx] domain
p.np=pde.grid.nPoints; p.nu=p.np; p.sol.xi=1/(p.nu); 
%[po,t,e]=getpte(p);
%p.mesh.bp=po; p.mesh.bt=t; p.mesh.be=e; % background mesh (for mesh adaption) 
h=par(1)*ones(p.np,1);
p.vol=2*lx;   % Total volume
par(2)=par(1)*(p.vol)/2;
p.u=h; p.u=[p.u; par']; % initial guess (homogeneous) 

p.nc.nsteps=20; p.sw.foldcheck=1; p.plot.auxdict={'lambda','c0','mu0'}; 
p.plot.pstyle=1; p.usrlam=[0 0.5 1]; p.nc.nsteps=100; p.sw.jac=1; 
p.sw.bifcheck=2;
p.fuha.outfu = @bra;      % Measured observables along the branch
p.spcontsw=0;   % Compute numerically spjac.m


% % Compute first momenta c1 in the initial state
% fem=p.pdeo.fem; gr=p.pdeo.grid; 
% x=getpte(p); x=x';
% [~,M,~]=fem.assema(gr,1,1,1);
% c1=sum(M*p.u(1:p.nu))/p.vol;
% par = p.u(p.nu+1:end);
% par(2)=c1;
% p.u=[p.u(1:p.nu); par];