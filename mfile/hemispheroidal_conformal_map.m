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



function map = hemispheroidal_conformal_map(v,f,c)
% A fast method for computing hemispheroidal conformal map of an open surface
%
% Input:
% v: nv x 3 vertex coordinates of a simply-connected triangle mesh
% f: nf x 3 triangulations of a simply-connected triangle mesh
% c: the z-radius of the target hemispheroid (the x,y radii are fixed as 1)
%
% Output:
% map: nv x 3 vertex coordinates of the hemispheroidal cap conformal parameterization
%
% Remark:
% - The input surface should be aligned with the xyz axis beforehand,
%   otherwise the result may be affected
%
% Written by Gary Pui-Tung Choi, 2023

nv = length(v);

% extract mesh boundary
bd = meshboundaries(f);
if length(bd)~=1
    % if the number of boundaries is not 1, the surface is not simply connected
    error('The input mesh is not a simply connected open surface!');
else
    bdy_v = bd{1};
end

%% Disk conformal map for the input surface

% Disk conformal map (using Choi and Lui, J. Sci. Comput. 2015)
disk = disk_conformal_map(v,f,bdy_v);

% (optional) double check the orientation of the disk conformal map 
e1 = [disk(f(:,2),:) - disk(f(:,1),:), disk(f(:,1),1).*0];
e2 = [disk(f(:,3),:) - disk(f(:,1),:), disk(f(:,1),1).*0];
cross12 = cross(e1,e2);
if sum(cross12(:,3) > 0) < sum(cross12(:,3) < 0)
    % wrong orientation in flattening
    disk = [disk(:,2), disk(:,1)];
    disp('Corrected orientation.');
end

% Search for an optimal Mobius transformation to further reduce the area distortion 
% (extending the method in Choi et al., SIAM J. Imaging Sci. 2020)

% Compute the area with normalization
area_v = face_area(f,v); 
area_v = area_v/sum(area_v);

z = complex(disk(:,1),disk(:,2));

% Fixing the rotation arbitrarily
[~,id_right] = max(v(:,1));
z = z*exp(-1i*angle(z(id_right)));

% Function for calculating the area after the Mobius transformation 
area_map = @(x) face_area(f,spheroidal_projection(...
    (z-x(1)*exp(1i*x(2)))./(1-conj(x(1)*exp(1i*x(2)))*z),1,c))/...
    sum(face_area(f,spheroidal_projection(...
    (z-x(1)*exp(1i*x(2)))./(1-conj(x(1)*exp(1i*x(2)))*z),1,c)));

% objective function: 
d_area = @(r) finitemean(abs(log(area_map(r)./area_v)).^2);

% Optimization setup
x0 = [0,0]; % initial guess, try something different if the result is not good
lb = [0,-pi]; % lower bound for the parameters
ub = [1,pi]; % upper bound for the parameters
options = optimoptions('fmincon','Display','off');

% Optimization
x = fmincon(d_area,x0,[],[],[],[],lb,ub,[],options);

% obtain the conformal parameterization with area distortion corrected
fz = (z-x(1)*exp(1i*x(2)))./(1-conj(x(1)*exp(1i*x(2)))*z);
disk_mobius = [real(fz),imag(fz)];

E_mobius = spheroidal_projection(fz,1,c);

%% Quasi-conformal composition

% compute Beltrami coefficient
mu = beltrami_coefficient(disk_mobius,f,E_mobius);

% Construct a QC map with the same Beltrami coefficient
bd = meshboundaries(f);
disk_map = linear_beltrami_solver(disk_mobius,f,mu,...
    bd{1},disk_mobius(bd{1},:));

% use disk_map to E_mobius to do conformal interpolation
Fx = scatteredInterpolant(disk_map(:,1),disk_map(:,2),E_mobius(:,1));
Fy = scatteredInterpolant(disk_map(:,1),disk_map(:,2),E_mobius(:,2));
Fz = scatteredInterpolant(disk_map(:,1),disk_map(:,2),E_mobius(:,3));

% Obtain the inverse mapping result
map = zeros(nv,3);
map(:,1) = Fx(disk_mobius(:,1),disk_mobius(:,2));
map(:,2) = Fy(disk_mobius(:,1),disk_mobius(:,2));
map(:,3) = Fz(disk_mobius(:,1),disk_mobius(:,2));

end


function m = finitemean(A)
    % for avoiding the Inf values caused by division by a very small area
    m = mean(A(isfinite(A)));
end
