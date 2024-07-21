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


% Demo for the proposed hemispheroidal parameterization methods:
% - Tutte: hemispheroidal_tutte_map.m
% - Conformal: hemispheroidal_conformal_map.m
% - Area-preserving: hemispheroidal_area_preserving_map.m
% - Balanced: hemispheroidal_balanced_map.m

addpath('input_geom');
addpath('mfile')

%% Load a mesh

fname = 'example_face.mat';
% fname = 'example_sophie.mat';
% fname = 'example_julius.mat';
% fname = 'example_chinese_lion.mat';
% fname = 'example_matterhorn.mat';
% fname = 'example_brain.mat';
% fname = 'example_bunny.mat';
load(fname);

% Optional: clean the mesh and remove all valence 2 vertices (vertices 
% connected to only 1 triangle) at the boundary
% [v,f] = clean_mesh(v,f,1); 

plot_mesh(v,f); 
view([20 20])

%% Determine a suitable hemispheroidal shape

% Transfer data to cononical coordinates (coordinates shift and rotation)
% Surface registration step:
[new_v, centroid, R, normal_sign] = register_surface(v, f, false);
% To restore the rotation of the coordinates use: (R * new_v)'

if normal_sign < 0
    new_v(:, 2:3) = new_v(:, 2:3) * normal_sign;
end

% Fit hemispheroidal shape over the mesh

% Bounding box after alignment
BB_d = max(new_v) - min(new_v);
BB_x = BB_d(1); 
BB_y = BB_d(2); 
BB_z = BB_d(3);
BB_xy = mean([BB_x, BB_y]);
BB_z_normalized = abs(BB_z) / BB_xy;
if BB_z_normalized >= 1
    hemispheroid_type = "prolate";
else
    hemispheroid_type = "oblate";
end
 

% Compute the spheroidal coordinates parameters (with a unit disk base)
foci = sqrt(abs(1 - BB_z_normalized^2)); % Focal distance
if hemispheroid_type == "prolate"
    zeta = asinh(1 / foci);
    aa = sinh(zeta) * foci; % "a" axis size (re-calculate)
    cc = cosh(zeta) * foci; % "c" axis size (re-calculate)

elseif hemispheroid_type == "oblate"
    zeta = acosh(1 / foci);
    aa = cosh(zeta) * foci; % "a" axis size (re-calculate)
    cc = sinh(zeta) * foci; % "c" axis size (re-calculate)

end
disp(strcat("The surface is ", hemispheroid_type, " with a = ",...
    num2str(aa), " and c = ", num2str(cc)))

%% Hemispheroidal tutte map
map_T = hemispheroidal_tutte_map(new_v,f,cc);

plot_mesh(map_T,f); 
view([20 20])

% Evaluate the angle and area distortion 
distortion_angle = angle_distortion(v,f,map_T);
distortion_area = area_distortion(v,f,map_T);

fprintf('Mean(angle distortion) = %.2f\n',mean(abs(distortion_angle)));
fprintf('SD(angle distortion) = %.2f\n',std(abs(distortion_angle)));
fprintf('Mean(area distortion) = %.2f\n',mean(abs(distortion_area)));
fprintf('SD(area distortion) = %.2f\n',std(abs(distortion_area)));


%% Hemispheroidal conformal map
map_C = hemispheroidal_conformal_map(new_v,f,cc);

plot_mesh(map_C,f); 
view([20 20])

% Evaluate the angle and area distortion 
distortion_angle = angle_distortion(v,f,map_C);
distortion_area = area_distortion(v,f,map_C);

fprintf('Mean(angle distortion) = %.2f\n',mean(abs(distortion_angle)));
fprintf('SD(angle distortion) = %.2f\n',std(abs(distortion_angle)));
fprintf('Mean(area distortion) = %.2f\n',mean(abs(distortion_area)));
fprintf('SD(area distortion) = %.2f\n',std(abs(distortion_area)));


%% Hemispheroidal area-preserving map
map_A = hemispheroidal_area_preserving_map(new_v,f,cc);

plot_mesh(map_A,f); 
view([20 20])

% Evaluate the angle and area distortion 
distortion_angle = angle_distortion(v,f,map_A);
distortion_area = area_distortion(v,f,map_A);

fprintf('Mean(angle distortion) = %.2f\n',mean(abs(distortion_angle)));
fprintf('SD(angle distortion) = %.2f\n',std(abs(distortion_angle)));
fprintf('Mean(area distortion) = %.2f\n',mean(abs(distortion_area)));
fprintf('SD(area distortion) = %.2f\n',std(abs(distortion_area)));


%% Hemispheroidal balanced map
alpha = 0.2; % Tutte weight
beta  = 0.3;  % conformal weight
gamma = 0.5; % area-preserving weight
map_B = hemispheroidal_balanced_map(v,f,cc,alpha,beta,gamma);

plot_mesh(map_B,f); 
view([20 20])

% Evaluate the angle and area distortion 
distortion_angle = angle_distortion(v,f,map_B);
distortion_area = area_distortion(v,f,map_B);

fprintf('Mean(angle distortion) = %.2f\n',mean(abs(distortion_angle)));
fprintf('SD(angle distortion) = %.2f\n',std(abs(distortion_angle)));
fprintf('Mean(area distortion) = %.2f\n',mean(abs(distortion_area)));
fprintf('SD(area distortion) = %.2f\n',std(abs(distortion_area)));

