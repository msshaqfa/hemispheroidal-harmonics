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


function map = hemispheroidal_tutte_map(v,f,c)
% A fast method for computing hemispheroidal Tutte map of a simply connected open surface
%
% Input:
% v: nv x 3 vertex coordinates of a simply connected triangle mesh
% f: nf x 3 triangulations of a simply connected triangle mesh
% c: the z-radius of the target hemispheroid (the x,y radii are fixed as 1)
%
% Output:
% map: nv x 3 vertex coordinates of the hemispheroidal Tutte map
%
% Remark:
% - The input surface should be aligned with the xyz axis beforehand,
%   otherwise the result may be affected
%
% Written by Gary Pui-Tung Choi, 2024

% extract mesh boundary
bd = meshboundaries(f);
if length(bd)~=1
    % if the number of boundaries is not 1, the surface is not simply connected
    error('The input mesh is not a simply connected open surface!');
else
    bdy_v = bd{1};
end

%% Disk Tutte map with arclength parameterization boundary constraint
bdy_length = sqrt((v(bdy_v,1) - v(bdy_v([2:end,1]),1)).^2 + ...
            (v(bdy_v,2) - v(bdy_v([2:end,1]),2)).^2 + ...
            (v(bdy_v,3) - v(bdy_v([2:end,1]),3)).^2);
partial_edge_sum = zeros(length(bdy_v),1);
for i = 2:length(bdy_length)
    for j = 1:i-1
    partial_edge_sum(i) = partial_edge_sum(i) + bdy_length(j);
    end
end
theta = 2*pi.*partial_edge_sum/sum(bdy_length)+2*pi*0.01;
bdy = exp(theta*1i)';
disk = tutte_map(v,f,bdy_v,bdy); 

% double check the orientation of the disk map 
e1 = [disk(f(:,2),:) - disk(f(:,1),:), disk(f(:,1),1).*0];
e2 = [disk(f(:,3),:) - disk(f(:,1),:), disk(f(:,1),1).*0];
cross12 = cross(e1,e2);
if sum(cross12(:,3) > 0) < sum(cross12(:,3) < 0)
    % wrong orientation in flattening
    disk = [disk(:,2), disk(:,1)];
    disp('Corrected orientation.');
end

%% Obtain the final hemispheroidal parameterization
map = spheroidal_projection(disk,1,c); 


end
