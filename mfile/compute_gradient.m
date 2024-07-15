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
% Authors of this file: Gary Choi @ 2024



function grad = compute_gradient(v,f,g)
% compute the gradient of each face
% input:
% v: k x 3 % vertex set
% f: m x 3 face set
% g: k x 1 a function defined on each vertex
% output: grad: m x 3
%
% If you use this code in your own work, please cite the following paper:
% [1] G. P. T. Choi and C. H. Rycroft, 
%     "Density-equalizing maps for simply connected open surfaces."
%     SIAM Journal on Imaging Sciences, 11(2), pp. 1134-1178, 2018.
% 
% Copyright (c) 2017-2018, Gary Pui-Tung Choi
% https://scholar.harvard.edu/choi

if size(v,2) ~= 3
    v = [v, v(:,1).*0];
end

% compute edges
e1 = v(f(:,3),:) - v(f(:,2),:);
e2 = v(f(:,1),:) - v(f(:,3),:);
e3 = v(f(:,2),:) - v(f(:,1),:);

% compute area
cross12 = cross(e1,e2);
area = abs(1/2*(cross12(:,1).^2+cross12(:,2).^2+cross12(:,3).^2).^(1/2));
N = [(1./(2*area)).*cross12(:,1), (1./(2*area)).*cross12(:,2), (1./(2*area)).*cross12(:,3)];

% compute gradient
temp = [g(f(:,1)).*e1(:,1),g(f(:,1)).*e1(:,2),g(f(:,1)).*e1(:,3)] + ...
    [g(f(:,2)).*e2(:,1),g(f(:,2)).*e2(:,2),g(f(:,2)).*e2(:,3)] + ...
    [g(f(:,3)).*e3(:,1),g(f(:,3)).*e3(:,2),g(f(:,3)).*e3(:,3)];
grad = cross(N,temp);
grad = [1./(2*area).*grad(:,1),1./(2*area).*grad(:,2)];

