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


% Demo for the proposed hemispheroidal harmonic decomposition framework

addpath('input_geom/')
addpath('mfile')
addpath('distmesh/')

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

%% Parameter setup

% Mapping type (uncomment one only)
% mapping_type = "tutte";
% mapping_type = "conformal";
mapping_type = "area_preserving";
% mapping_type = "balanced"; alpha = 0.2; beta = 0.3; gamma = 0.5;

% Harmonic reconstruction setting
max_n = 50; % Max Analysis degree
rec_max_n = max_n; % Max reconstruction degree
edge_length = 0.0250;% Reconstruction domain resolution


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

%% Hemispheroidal parameterization

switch mapping_type
    case "tutte"
        map = hemispheroidal_tutte_map(new_v, f, cc);
    case "conformal"
        map = hemispheroidal_conformal_map(new_v, f, cc);
    case "area_preserving"
        map = hemispheroidal_area_preserving_map(new_v, f, cc);
    case "balanced"
        map = hemispheroidal_balanced_map(new_v, f, cc, alpha, beta, gamma);
end

plot_mesh(map,f);
view([20 20])

%% Harmonic expansion (Solver)
[thetas, phis] = cart2spheroid(map, foci, zeta, hemispheroid_type);
D_mat = hemispheroidal_harmonic_basis(max_n, thetas, phis, hemispheroid_type);
qm_k = D_mat\new_v;

%% Harmonic reconstruction
epsillon = pi/50;  % To account for numerical errors.

theta_c = pi /2 - epsillon;
reconstruction_mesh_option = 1;

% Construct reconstruction domain
[rec_v, rec_f] = uniform_spherical_cap_grid(theta_c, edge_length, 100, reconstruction_mesh_option);
rec_v(:, 3) = rec_v(:, 3) * cc; % This will slightly deform the areas/angles of the triangles

[rec_thetas, rec_phis] = cart2spheroid(rec_v, foci, zeta, hemispheroid_type);

rec_D_mat = hemispheroidal_harmonic_basis(rec_max_n, rec_thetas, ...
    rec_phis, hemispheroid_type);

rec_v = real(rec_D_mat * qm_k(1:(rec_max_n+1)^2, :));

rec_v_old = rec_v;

% To restore the surface registration
rec_v(:, 2:3) = normal_sign .* rec_v(:, 2:3);
rec_v = (R * rec_v')';
rec_v = rec_v + centroid;

plot_mesh(rec_v, rec_f);
view([20 20])

