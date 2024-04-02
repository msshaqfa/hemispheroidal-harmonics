function D_mat = hemispheroidal_harmonic_basis_math(max_n, thetas, phis, hemispheroid_type)
% Imeplementation: Mahmoud Shaqfa.

if hemispheroid_type == "oblate"
    xi = 2. * sin(thetas) - 1;
elseif hemispheroid_type == "prolate"
    xi = 1 - 2. * cos(thetas);
end

D_mat = zeros([size(thetas), (max_n+1)^2]);

Legendre_table = legendre(max_n, xi)';

for n = 0:max_n
%   For positive orders
    for m = 0:n
        % Neumann BCs
        Normalization = sqrt((2*n +1) * factorial(n-m) / (4*pi * factorial(n+m)));
        D_mat(:, n^2 + n + m + 1) = Normalization .* Legendre_table(:, m+1) .* exp(1i*m.*phis);
    end
%   For redundunt negative orders
    for m = -n:-1
        D_mat(:, n^2 + n + m + 1) = conj(D_mat(:, n^2 + n + abs(m) + 1)) .* (-1)^abs(m);
    end
    fprintf("\nCompleted: %3.2f %% for n = %3.0f", ((n+1)^2/(max_n+1)^2*100), n)
end