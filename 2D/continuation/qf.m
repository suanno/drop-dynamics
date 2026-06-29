function q = qf(p,u)
%% q: array of aux. conditions as defined below

%{

Here, q has only one entry for mean concentration conservation
c0: mean concentration (par# 2)

%}

%%
par = u(p.nu+1:end);
u = u(1:p.nu);
x=getpte(p); x=x';
c1 = par(2);

% Mass matrix (needed for (riemann) integrals)
M=p.mat.M;

% Integral constraint
q = sum(M*(x.*u))/p.vol - c1; 
end