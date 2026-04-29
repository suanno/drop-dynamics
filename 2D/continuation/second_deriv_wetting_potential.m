function out = second_deriv_wetting_potential(h, h_a)
    % Derivative of the potential, WITHOUT the minus sign in front
    out = 6*h_a.^3*h.^(-7) - 3*h.^(-4);
end