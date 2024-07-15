%*************************************************************************%
% All rights reserved (C) to the authors: Mahmoud SHAQFA, Gary CHOI       %
%                                                                         %
%                                                                         %
% M. Shaqfa Contact:                                                      %
% Department of mechanical engineering,                                   %
% Massachusetts Institute of Technology (MIT)                             %
% Cambridge, MA, USA                                                      %
%               Email: mshaqfa@mit.edu                                    %
% G. Choi Contact:                                                        %
% Department of Mathematics, The Chinese University of Hong Kong,         %
% Hong Kong                                                               %
%               Email: ptchoi@cuhk.edu.hk                                 %
%                                                                         %
%*************************************************************************%
% This code includes implementations for:                                 %
%				- Hemispheroidal harmonics analysis (HSOH)                %
%				- Hemispheroidal parameterization                         %
%                 -- Area-preserving mapping                              %
%                 -- Conformal mapping                                    %
%                 -- Tutte mapping                                        %
%                 -- Balanced mapping                                     %
% This code is part of the paper: "Hemispheroidal parameterization and    %
% harmonic decomposition of simply connected open surfaces"               %
%                                                                         %
%*************************************************************************%
% This library is free software; you can redistribute it and/or modify	  %
% it under the terms of the GNU Lesser General Public License as published%
% by the Free Software Foundation; either version 2.1 of the License, or  %
% (at your option) any later version.                                     %
%                                                                         %
% This library is distributed in the hope that it will be useful,         %
% but WITHOUT ANY WARRANTY; without even the implied warranty of          %
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                    %
% See the GNU Lesser General Public License for more details.        	  %
% You should have received a copy of the GNU Lesser General Public License%
% along with this library; if not, write to the Free Software Foundation, %
% Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA       %
%*************************************************************************%
% Authors of this file: Mahmoud S. Shaqfa @ 2024



function D_mat = hemispheroidal_harmonic_basis(max_n, thetas, phis, hemispheroid_type, varargin)
% Imeplementation: Mahmoud Shaqfa.

if nargin == 4
    printed = true;
else
    printed = varargin{1};
end

if hemispheroid_type == "oblate"
    xi = 2. * sin(thetas) - 1;
elseif hemispheroid_type == "prolate"
    xi = 1 - 2. * cos(thetas);
end

% To account for numerical round up error
% (e.g., -1.00001 or 1.00002 are not acceptable as \xi \in [-1, 1])
epsillon = 10e-3;
xi(xi >= 1) = 1 - epsillon;
xi(xi <= -1) = -1 + epsillon;

D_mat = zeros(size(thetas, 1), (max_n+1)^2);

for n = 0:max_n
    Legendre_table = legendre(n, xi)';
%   For positive orders
    for m = 0:n
        % Neumann BCs
        Normalization = sqrt((2*n +1) * factorial(n-m) / (4 * pi * factorial(n+m)));
        D_mat(:, n^2 + n + m + 1) = Normalization .* Legendre_table(:, m+1) .* exp(1i .* m.* phis);
    end
%   For redundunt negative orders
    for m = -n:-1
        D_mat(:, n^2 + n + m + 1) = conj(D_mat(:, n^2 + n + abs(m) + 1)) .* (-1) .^ abs(m);
    end
    
    if printed
        fprintf("\nCompleted: %3.2f %% for n = %3.0f", ((n+1)^2/(max_n+1)^2*100), n)
    end
end
