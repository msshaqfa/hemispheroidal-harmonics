function D_mat = disk_harmonic_math_reconstruction(K, eigen_table, reconstruction_resolution, qm_k)
% Imeplementation: Mahmoud Shaqfa, Gary Choi, Katrin Beyer.
% Bessel-Fourier basis functions for the disk harmonics (\rho, \phi).
% It returns a table for positive and negative orders for a certain index K
% The number of the columns in the matrix is 2K+1 for (-m:1:m) orders
% Note: rhos and phis matrix msut be identical in size
% All the angles here must be in radians.
% For our applications the following variables shall be fixed:
% rhos and phis must be square matrix or column matrix (row matrix will give you an error!)
% BC = "Neumann";
[phis, rhos] = meshgrid(linspace(0, 2*pi, reconstruction_resolution),...
    linspace(0, 1, reconstruction_resolution));
sz = size(rhos);
sz(3) = (K+1)^2;
D_mat = zeros(sz);
temp_sz = size(rhos); temp_sz(3) = 3;
reconstruction = zeros(temp_sz);
for n = 0:K
    for m = 0:n
        lmn = eigen_table(n+1, abs(m)+1);
        Normalization = 0.5 * (1 - (abs(m)/lmn)^2) * besselj(abs(m), lmn)^2;
        Normalization = 1/sqrt(Normalization);          % Bessel normalization (for Neumann BC only)
        Normalization = Normalization * (1/sqrt(2*pi)); % Fourier normalization
        D_mat(:, :, n^2 + n + m + 1) = Normalization .* besselj(abs(m), lmn*rhos) .* exp(1i*abs(m).*phis);
    end
    for m = -n:-1
        D_mat(:, :, n^2 + n + m + 1) = conj(D_mat(:, :, n^2 + n + abs(m) + 1)) .* (-1)^m;
    end
end
for n = 0:K
    for ii = n^2 + 1:n^2 + 2*n + 1
        reconstruction(:, :, 1) = reconstruction(:, :, 1) + D_mat(:, :, ii) .* qm_k(ii, 1);
        reconstruction(:, :, 2) = reconstruction(:, :, 2) + D_mat(:, :, ii) .* qm_k(ii, 2);
        reconstruction(:, :, 3) = reconstruction(:, :, 3) + D_mat(:, :, ii) .* qm_k(ii, 3);
    end
    fprintf("\nCompleted: %3.2f %% for n = %3.0f", ((n+1)^2/(K+1)^2*100), n)
end