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
