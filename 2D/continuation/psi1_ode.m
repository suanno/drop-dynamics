function dpsidr = psi1_ode(r, psi, Qfun, hfun, Qrfun, hrfun)
    % Inputs are the functions Q(r) and h(r) (and first derivatives) found through interpolations
    dpsidr = [psi(2); -(1/r)*psi(2)*(1-(r/Qfun(r))*Qrfun(r))+psi(2)/r-hrfun(r)];
end