function q = qf(p,u)
%% q: array of aux. conditions as defined below

%{

Here, q has only one entry for mean concentration conservation
c0: mean concentration (par# 2)

%}

%%
par = u(p.nu+1:end);
up = u(1:p.nu);
c0 = par(2);

% Mass matrix (needed for (riemann) integrals)
M0 = p.mat.M0;
u=p.mat.fill*up;

% Integral constraint: q = 1/vol*integral(u) - c0 := 0
q1 = sum(M0*u)/p.vol - c0; %The domain size is par(5)

% Phase condition (translational invariance)
%uold=p.u(1:p.nu); u0x=p.mat.Kx*uold; q2=u0x'*u(1:p.nu);

% Center of mass condition
x=getpte(p); x=x'; % extract the point coordinates from p
q2 = sum(M0*(u.*x(1:end)))-p.xcm;

% To activate
if p.reducedmass > 0
    % Reduced mass (mass ABOVE the wetting layer)
    q1 = sum(M0*(u-min(u)))/p.vol - p.reducedmass;
end

switch p.nc.nq
    case 1; q=q1;
    case 2; q=[q1;q2];
end