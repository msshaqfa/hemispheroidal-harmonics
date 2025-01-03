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


function [new_v, centroid, V, normal_sign] = register_surface(v, f, plot_figs)
% Written by: Mahmoud Shaqfa
% Usage:
% [new_V, centroid, r] = register_surface(v, verts_limit)
%
% Input:
% v: nv x 3 vertex coordinates of a simply-connected triangle mesh
% verts_limit: the number of verts used in the registration process (randomly sampled)
% plot_figs: true or false for plotting the boundaries with the fitted
% plane.

% Extract mesh boundaries:
[~, outer] = meshboundaries(f);
bds_v = v(outer, :);

% First step: Center the surface data around the origin (considers only the
% center of mass for the boundaries in this work)
centroid = mean(bds_v); % Centroid of boundaries
new_v = v - centroid;
bds_v =  bds_v - centroid;

bds_v_cpy = bds_v;


% Second step: Apply SVD for the surface data (only the exterior boundary)
C_plane = fit_plane(bds_v);

norm_plane = sqrt(sum(C_plane.^2));


% Compute angles for tests
% l = C_plane(1) / norm_plane;
% m = C_plane(2) / norm_plane;
n = C_plane(3) / norm_plane;
normal_sign = sign(-n);

% tx = acosd(l);
% ty = acosd(m);
% tz = acosd(n);

% res = 10; % Just any small number > 2
res = 30; % Just any small number > 2
plane_x = linspace(min(bds_v(:, 1)), max(bds_v(:, 1)), res);
plane_y = linspace(min(bds_v(:, 2)), max(bds_v(:, 2)), res);
[plane_x, plane_y] = meshgrid(plane_x, plane_y);
plane_z = C_plane(1) .* plane_x + C_plane(2) .* plane_y + C_plane(3); % The plane equation

plane_x_vec = reshape(plane_x, 1, []);
plane_y_vec = reshape(plane_y, 1, []);
plane_z_vec = reshape(plane_z, 1, []);
plane_scatter = [plane_x_vec; plane_y_vec; plane_z_vec]';

% [U, Sigma, V] = svd(bds_v);
[U, Sigma, V] = svd(plane_scatter);

% Check the normal (flip if necessary)
if normal_sign < 0
    disp("Flip the normal!")
    rotz = @(t) [cos(t) -sin(t) 0 ; sin(t) cos(t) 0 ; 0 0 1];
    V = V * rotz(pi()/5);
end

new_v = (V' * new_v')';
new_bds_v = (V' * bds_v')';


if plot_figs
    scatter3(new_bds_v(:, 1), new_bds_v(:, 2), new_bds_v(:, 3), 'kx', 'DisplayName', 'New boundary')
    % plot3(new_bds_v(:, 1), new_bds_v(:, 2), new_bds_v(:, 3), 'red', 'LineWidth', 4)
    
    hold on
    scatter3(bds_v_cpy(:, 1), bds_v_cpy(:, 2), bds_v_cpy(:, 3), 'ro', 'DisplayName', 'Old boundary')
    % plot3(bds_v_cpy(:, 1), bds_v_cpy(:, 2), bds_v_cpy(:, 3), 'blue', 'LineWidth', 4)
    
    hold on
    surf(plane_x, plane_y, plane_z, 'DisplayName', 'Plane fit', 'AlphaData', 0.1,...
        'FaceAlpha', 0.5, 'EdgeColor', 'none')
    hold on
    C_plane2 = fit_plane(new_bds_v); % second iteration for the new boundaries
    plane_z2 = C_plane2(1) .* plane_x + C_plane2(2) .* plane_y + C_plane2(3); % The plane equation
    surf(plane_x, plane_y, plane_z2, 'DisplayName', 'Plane fit (final)', 'AlphaData', 0.1,...
        'FaceAlpha', 0.2, 'EdgeColor', 'none', 'FaceColor','b')
    hold off
    legend
    axis equal tight off
    box on;
    title("Boundary fitting for surface registratrion")
    xlabel('x')
    ylabel('y')
    zlabel('z')
end

disp(strcat("The normal sign is: ", num2str(normal_sign)))
end


function C = fit_plane(v)
% A modified version by: Mahmoud Shaqfa
% This is to fit the implicit equation of a first order 2D plane.
% Original Source: Val Schmidt (2023). planefit (https://www.mathworks.com/matlabcentral/fileexchange/36353-planefit), 
% MATLAB Central File Exchange.
% Retrieved November 27, 2023.

xx = v(:, 1);
yy = v(:, 2);
zz = v(:, 3);
O = ones(length(xx), 1);
C = [xx yy O]\zz; % Least squares
end
