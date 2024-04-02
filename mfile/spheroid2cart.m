function cart_coords = spheroid2cart(theta, phi, foci, zeta, spheroid_type)

if spheroid_type == "prolate"
    x = foci .* sinh(zeta) .* sin(theta) .* cos(phi);
    y = foci .* sinh(zeta) .* sin(theta) .* sin(phi);
    z = foci .* cosh(zeta) .* cos(theta);

elseif spheroid_type == "oblate"
    x = foci .* cosh(zeta) .* cos(theta) .* cos(phi);
    y = foci .* cosh(zeta) .* cos(theta) .* sin(phi);
    z = foci .* sinh(zeta) .* sin(theta);
end

cart_coords = [x, y, z];
end