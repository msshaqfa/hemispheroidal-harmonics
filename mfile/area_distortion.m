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



function distortion = area_distortion(v,f,map)

% Calculate and visualize the area distortion log(area_map/area_v)
% 
% Input:
% v: nv x 3 vertex coordinates of a genus-0 triangle mesh
% f: nf x 3 triangulations of a genus-0 triangle mesh
% map: nv x 2 or 3 vertex coordinates of the mapping result
%
% Output:
% distortion: 3*nf x 1 area differences
% 
% If you use this code in your own work, please cite the following paper:
% [1] G. P. T. Choi, Y. Leung-Liu, X. Gu, and L. M. Lui, 
%     "Parallelizable global conformal parameterization of simply-connected surfaces via partial welding."
%     SIAM Journal on Imaging Sciences, 2020.
%
% Copyright (c) 2018-2020, Gary Pui-Tung Choi
% https://scholar.harvard.edu/choi

nv = length(v);
nv2 = length(map);

if nv ~= nv2
    error('Error: The two meshes are of different size.');
end

if size(v,2) == 1
    v = [real(v),imag(v),zeros(length(v),1)];
elseif size(v,2) == 2
    v = [v,zeros(length(v),1)];
end

if size(map,2) == 1
    map = [real(map),imag(map),zeros(length(map),1)];
elseif size(map,2) == 2
    map = [map,zeros(length(map),1)];
end

% calculate area of v
area_v = face_area(f,v);
% calculate area of map
area_map = face_area(f,map);

% normalize the total area
v = v*sqrt(sum(area_map)/sum(area_v));
area_v = face_area(f,v);

% calculate the area ratio
distortion = log(area_map./area_v);
