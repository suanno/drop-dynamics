function out = wetting_potential(H) 
    % H=h/h_a is the adimensional height
    out = H.^(-5)/5 - H.^(-2)/2;
end
