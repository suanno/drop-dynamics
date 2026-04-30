function out = bra(p, u)
h = u(1:p.nu); % Exclude the parameters of the pde
par = u(p.nu+1:end); % Parameters
x=getpte(p); x=x';
e1 = ones(1,length(x))';

ha = par(1);
mass = par(2);  % Integral of h1 (without the hat!)
hout = min(h);  %Measured hout
hhat = h-hout*e1;
what = wetting_potential(h,ha)-wetting_potential(hout,ha)*e1;
Qin = h.^3;


%M0 = p.mat.M0;

% Integral constraint: q = 1/vol*integral(u) - c0 := 0
%a2 = sum(M*hhat);


a1 = trapz(x, what);
a2 = trapz(x, hhat);
a3 = trapz(x, hhat./Qin);

out = [ha; mass; hout; a1; a2; a3];