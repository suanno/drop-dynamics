function out = bra(p, u)
h = u(1:p.nu); % Exclude the parameters of the pde
par = u(p.nu+1:end); % Parameters
x=getpte(p); x=x';
e1 = ones(1,length(x))';

mass = par(1);  % Integral of h1 (without the hat!)
hout = min(h);  %Measured hout
hhat = h-hout*e1;
what = wetting_potential(h)-wetting_potential(hout)*e1;
Qin = h.^3;%Should be divided by 3??


%M0 = p.mat.M0;

% Integral constraint: q = 1/vol*integral(u) - c0 := 0
%a2 = sum(M*hhat);

% Copy from 1D!!!
a1 = trapz(x, what);
a2 = trapz(x, hhat);
a3 = trapz(x, hhat./Qin);

out = [mass; hout; a1; a2; a3];