function out = bra(p, u)
H = p.mat.fill*u(1:p.nu); % Exclude the parameters of the pde
par = u(p.nu+1:end); % Parameters
x=getpte(p); x=x';
e1 = ones(1,length(x))';

c0 = par(1);  % Integral of h1 (without the hat!)

%hout from hmax
Hmax = max(H);
Hout = 1+(wetting_potential(Hmax)-wetting_potential(1))/((Hmax-1)*second_deriv_wetting_potential(1));
%hout = min(h);  %Measured hout
Hhat = H-Hout*e1;
omegahat = wetting_potential(H)-wetting_potential(Hout)*e1;

%M0 = p.mat.M0;

% Integral constraint: q = 1/vol*integral(u) - c0 := 0
%a2 = sum(M*hhat);

%dhomegaout = par(3);
Omegahat = sum(p.mat.M0*(omegahat));
m = sum(p.mat.M0*Hhat);
I = sum(p.mat.M0*(Hout^3.*Hhat./H.^3));
K = sum(p.mat.M0*(Hhat.^2./H.^3));

out = [c0; Hmax; Hout; m; Omegahat; I; K];