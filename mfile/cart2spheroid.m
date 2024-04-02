function [theta, phi] = cart2spheroid(v, foci, zeta, spheroid_type)

if spheroid_type == "prolate"
    xx = v(:, 1) ./ (foci .* sinh(zeta));
    yy = v(:, 2) ./ (foci .* sinh(zeta));
    zz = v(:, 3) ./ (foci .* cosh(zeta));
    
    theta = atan2(zz, (sqrt(xx .^ 2 + yy .^ 2))); 
    theta = pi/2 - theta;

elseif spheroid_type == "oblate"
    xx = v(:, 1) ./ (foci .* cosh(zeta));
    yy = v(:, 2) ./ (foci .* cosh(zeta));
    zz = v(:, 3) ./ (foci .* sinh(zeta));

    theta = atan2(zz, (sqrt(xx .^ 2 + yy .^ 2)));
end

phi = atan2(yy, xx);

end