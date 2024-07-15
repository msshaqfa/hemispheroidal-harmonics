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



function M = f2v_area(v,f)
% Face to vertex interpolation with area weighting.
%
% If you use this code in your own work, please cite the following paper:
% [1] G. P. T. Choi and C. H. Rycroft, 
%     "Density-equalizing maps for simply connected open surfaces."
%     SIAM Journal on Imaging Sciences, 11(2), pp. 1134-1178, 2018.
% 
% Copyright (c) 2017-2018, Gary Pui-Tung Choi
% https://scholar.harvard.edu/choi

nv = length(v);
nf = length(f);

if size(v,2) == 2
    v = [v,zeros(nv,1)];
end

% find area
area = face_area(f,v);

% create matrix
row = [f(:,3); f(:,1); f(:,2)];
col = [1:length(f), 1:length(f), 1:length(f)]';
val = [area; area; area];
M = sparse(row,col,val,nv,nf);
 
% normalize
vertex_area_sum = sum(M,2);
[Mrow,Mcol,Mval] = find(M);
M = sparse(Mrow,Mcol,Mval./vertex_area_sum(Mrow),nv,nf);
