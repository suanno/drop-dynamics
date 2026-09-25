function out = bra(p, u)
%H = p.mat.fill*u(1:p.nu); % Exclude the parameters of the pde
H=u(1:p.nu);
par = u(p.nu+1:end); % Parameters
x=getpte(p); x=x';
e1 = ones(1,length(x))';

c0 = par(1);  % Integral of h1 (without the hat!)

%hout from hmax
Hmax = max(H);
Hout_lim = min(H);
Hout = 1+(wetting_potential(Hmax)-wetting_potential(1))/((Hmax-1)*second_deriv_wetting_potential(1));
%hout = min(h);  %Measured hout
Hhat = H-Hout*e1;
omegahat = wetting_potential(H)-wetting_potential(Hout)*e1;

%M0 = p.mat.M0;

% Integral constraint: q = 1/vol*integral(u) - c0 := 0
%a2 = sum(M*hhat);

%dhomegaout = par(3);


% All integrals are multiplied by a factor 2, as with nbc we're solving half the
% droplet's shape (and so integrating half the space)
Omegahat = 2*sum(p.mat.M*(omegahat));
V = 2*sum(p.mat.M*Hhat);
I = 2*sum(p.mat.M*(Hout^3.*Hhat./H.^3));
K = 2*sum(p.mat.M*(Hhat.^2./H.^3));

out = [c0; Hmax; Hout_lim; Hout; V; Omegahat; I; K];