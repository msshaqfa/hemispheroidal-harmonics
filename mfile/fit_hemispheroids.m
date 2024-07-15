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



function fit_solution = fit_hemispheroids(v)
% This function to fit a hemispheroidal shape of the input surface

x_ = v(:, 1);
y_ = v(:, 2);
z_ = v(:, 3);

sz = size(x_);
A = zeros([sz(1), 3]);

A(:, 1) = x_.^ 2;
A(:, 2) = y_.^ 2;
A(:, 3) = z_.^ 2;
b = ones([size(x_), 1]);
fit_solution = lsqr(A, b);

c_ellips = [sqrt(1/fit_solution(1)), sqrt(1/fit_solution(2)), sqrt(1/fit_solution(3))];

disp("1st Stage of fitting ellispoid with a, b, and c:")
disp(c_ellips)


% Classify the hemispehroid
dist = zeros([3, 1]);
dist(1) = abs(c_ellips(1) - c_ellips(2));
dist(2) = abs(c_ellips(2) - c_ellips(3));
dist(3) = abs(c_ellips(1) - c_ellips(3));
min_dist = find(dist == min(dist));

disp(min_dist)

if min_dist == 1
    aa = mean([c_ellips(1), c_ellips(2)]);
    cc = c_ellips(3);
elseif min_dist == 2
    aa = mean([c_ellips(2), c_ellips(3)]);
    cc = c_ellips(1);
elseif min_dist == 3
    aa = mean([c_ellips(1), c_ellips(3)]);
    cc = c_ellips(2);
end

if aa <= cc
    hemispheroid_type = "prolate";
else
    hemispheroid_type = "oblate";
disp(strcat("The surface is ", hemispheroid_type, " with a size of a = ", num2str(aa), " and c = ", num2str(cc) ))

% 
% f = @(x,y,z) x.^2 + y.^2 - z.^2;
% interval = [-5 5 -5 5 0 5];
% fimplicit3(f)


end
