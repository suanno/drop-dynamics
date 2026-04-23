function out = bra(p, u)
h = p.mat.fill*u(1:p.nu); % Exclude the parameters of the pde
par = u(p.nu+1:end); % Parameters
x=getpte(p); x=x';
e1 = ones(1,length(x))';

ha = par(1);
c0 = par(2);  % Integral of h1 (without the hat!)

%hout from hmax
hmax = max(h);
hout = ha+(wetting_potential(hmax,ha)-wetting_potential(ha,ha))/((hmax-ha)*second_deriv_wetting_potential(ha,ha));
%hout = min(h);  %Measured hout
hhat = h-hout*e1;
what = wetting_potential(h,ha)-wetting_potential(hout,ha)*e1;

%M0 = p.mat.M0;

% Integral constraint: q = 1/vol*integral(u) - c0 := 0
%a2 = sum(M*hhat);

dWout = par(3);
Omega = sum(p.mat.M0*(what));
I = sum(p.mat.M0*(hout^3.*hhat./h.^3));
K = sum(p.mat.M0*(hhat.^2./h.^3));

out = [ha; c0; hmax; dWout; hout; Omega; I; K];