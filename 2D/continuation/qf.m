function q = qf(p,u)
%% q: array of aux. conditions as defined below

%{

Here, q has only one entry for mean concentration conservation
c0: mean concentration (par# 2)

%}

%%
par = u(p.nu+1:end);
u = u(1:p.nu);
c0 = par(2);

% Mass matrix (needed for (riemann) integrals)
M = p.mat.M;

% Integral constraint: q = 1/vol*integral(u) - c0 := 0
q = sum(M*u)/p.vol - c0; %The domain size is par(5)
end