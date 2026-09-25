function out = deriv_wetting_potential(H)
    out = -H.^(-6)+H.^(-3);
end