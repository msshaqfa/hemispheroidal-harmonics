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



function L = cotangent_laplacian(v,f)
% Compute the cotagent Laplacian of a mesh.
%
% If you use this code in your own work, please cite the following paper:
% [1] P. T. Choi and L. M. Lui, 
%     "Fast Disk Conformal Parameterization of Simply-Connected Open Surfaces."
%     Journal of Scientific Computing, 65(3), pp. 1065-1090, 2015.
%
% Copyright (c) 2014-2018, Gary Pui-Tung Choi
% https://scholar.harvard.edu/choi


nv = length(v);

f1 = f(:,1); f2 = f(:,2); f3 = f(:,3);

l1 = sqrt(sum((v(f2,:) - v(f3,:)).^2,2));
l2 = sqrt(sum((v(f3,:) - v(f1,:)).^2,2));
l3 = sqrt(sum((v(f1,:) - v(f2,:)).^2,2));

s = (l1 + l2 + l3)*0.5;
area = sqrt( s.*(s-l1).*(s-l2).*(s-l3));
 
cot12 = (l1.^2 + l2.^2 - l3.^2)./area/2;
cot23 = (l2.^2 + l3.^2 - l1.^2)./area/2; 
cot31 = (l1.^2 + l3.^2 - l2.^2)./area/2; 
diag1 = -cot12-cot31; diag2 = -cot12-cot23; diag3 = -cot31-cot23;

II = [f1; f2; f2; f3; f3; f1; f1; f2; f3];
JJ = [f2; f1; f3; f2; f1; f3; f1; f2; f3];
V = [cot12; cot12; cot23; cot23; cot31; cot31; diag1; diag2; diag3];
L = sparse(II,JJ,V,nv,nv);

end
