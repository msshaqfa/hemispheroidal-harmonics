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
% Authors of this file: Mahmoud S. Shaqfa and Gary Choi @ 2024

% Bijection between the unit disk and the Northern hemispheroid with radii a,a,c.
function P = spheroidal_projection(p,a,c)
    if size(p, 2) == 1
      p = [real(p), imag(p)];
    end
    u = p(:,1);
    v = p(:,2);
    if size(p,2) < 3
      z = 1 + u.^2 + v.^2;
      P = [2*a*u./z, 2*a*v./z, -c*(-1+u.^2+v.^2)./z];
      
      P(isnan(z)|(~isfinite(z)),1) = 0;
      P(isnan(z)|(~isfinite(z)),2) = 0;
      P(isnan(z)|(~isfinite(z)),3) = -c;
        
    else
      z = p(:,3);
      P = [(u/a)./(1+z/c), (v/a)./(1+z/c)];
      P(isnan(P)) = Inf;
    end
end
