function out = wetting_potential(h, h_a) 
    % h_a is the "heigh-scale"
    out = h_a^3/5*h.^(-5) - 1/2*h.^(-2);
end
