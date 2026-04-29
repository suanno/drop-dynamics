function out = deriv_wetting_potential(h, h_a)
    % Derivative of the potential, WITHOUT the minus sign in front
    out = - h_a.^3*h.^(-6) + h.^(-3);
end