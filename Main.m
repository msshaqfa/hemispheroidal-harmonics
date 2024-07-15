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


close all; clear; clc;

% addpath('code')
addpath('mfile')
addpath('stlTools')
addpath('distmesh/')
addpath('input_geom/')
% addpath('input_geom/others/')
addpath('rec_output/')

%% Load meshes and files

register = true;

 % For loading stl files:
% fname = "3D_face_refined.stl"; % (n=40 from the SCHA paper)
% fname = "Matterhorn_new_mode.stl"; % (n=60 from the DHA paper)
fname = "half_stone.stl";
% fname = "sophie.stl";
% fname = "bunny_open.stl";

[v, f, ~, ~] = stlRead(fname); % Uncomment only when loading STL files.

% For loading *.mat files
% fname = 'human_brain.mat';

% load(fname);  % Uncomment only when loading MAT files.


% fname = "bunny_open.off";
% [v,f] = read_off(fname);
% v=v'; f=f';

% Load obj files
% fname = "sophie.obj";
% obj_surf = readObj(fname);
% v = obj_surf.v;
% f = obj_surf.f;
% texture = obj_surf.vt;





% Mapping type (uncomment one only)
% mapping_type = "tutte";
% mapping_type = "conformal";
mapping_type = "area_preserving";
% mapping_type = "balanced";
alpha = 0.35;
beta = 0.0;
gamma = 0.65;


% Parameter setting
max_n = 10; % Max Analysis degree
rec_max_n = max_n; % Max reconstruction degree
edge_length = 0.0250;% Reconstruction domain resolution

plot_mesh(v,f);
view([-15 15]); 
axis on;
xlabel('x')
ylabel('y')
zlabel('z')




%% Transfer data to cononical coordinates (coordinates shift and rotation)
% Surface registration step:
[new_v, centroid, R, normal_sign] = register_surface(v, f, true);

if max(new_v(:, 3)) ~= max(abs(new_v(:, 3)))
    new_v(:, 3) = new_v(:, 3) * -1;
end

plot_mesh(new_v, f);

view([-15 15]); 
axis on;
xlabel('x')
ylabel('y')
zlabel('z')


%% Fit hemispheroidal shape over the mesh
% This part needs more work (Discuss with Gary the optimization)

% Bounding box after alignment
BB_d = max(new_v) - min(new_v);
BB_x = BB_d(1); BB_y = BB_d(2); BB_z = BB_d(3);

BB_xy = mean([BB_x, BB_y]);
BB_z_normalized = abs(BB_z) / BB_xy;

if BB_z_normalized >= 1
    hemispheroid_type = "prolate";
else
    hemispheroid_type = "oblate";
end

%% Compute the spheroidal coordinates parameters (with a unit disk base)
foci = sqrt(abs(1 - BB_z_normalized^2)); % Focal distance (of normalized base r = 1)
if hemispheroid_type == "prolate"
    zeta = asinh(1 / foci);
    aa = sinh(zeta) * foci; % "a" axis size (re-calculate)
    cc = cosh(zeta) * foci; % "c" axis size (re-calculate)

    % Plot implicit surface
    E_implicit = @(x, y, z) (x.^2 + y.^2) ./ (foci * sinh(zeta))^2 + z.^2 ./ (foci * cosh(zeta))^2 - 1;
    interval = [-1 1 -1 1 0 foci * cosh(zeta)];
    fimplicit3(E_implicit, interval,'EdgeColor','none','FaceAlpha',.5)
    title("Fitted prolate hemispheroid")
    axis equal

elseif hemispheroid_type == "oblate"
    zeta = acosh(1 / foci);
    aa = cosh(zeta) * foci; % "a" axis size (re-calculate)
    cc = sinh(zeta) * foci; % "c" axis size (re-calculate)

    % Plot implicit surface
    E_implicit = @(x, y, z) (x.^2 + y.^2) ./ (foci * cosh(zeta))^2 + z.^2 ./ (foci * sinh(zeta))^2 - 1;
    interval = [-1 1 -1 1 0 foci * sinh(zeta)];
    fimplicit3(E_implicit, interval,'EdgeColor','none','FaceAlpha',.5)
    view([-21.1000, 24.1389])
    axis equal off
end
disp(strcat("The surface is ", hemispheroid_type, " with a = ",...
    num2str(aa), " and c = ", num2str(cc)))
%% Ellipsoidal cap conformal map

% Hemispheroidal map 
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
fname_map_out = strcat('./rec_output/map_', mapping_type, '_', fname(1:end));
if endsWith(fname_map_out, '.mat', 'IgnoreCase', true)
    fname_map_out = strrep(fname_map_out, '.mat', '.stl');
end
stlWrite(fname_map_out, f, map);

% plot_mesh(map,f,mean_curv); view([-180 90]); 
plot_mesh(map,f); view([-180 90]); 

axis on;
xlabel('x')
ylabel('y')
zlabel('z')

% Evaluate the angle and area distortion 
distortion_angle = angle_distortion(v,f,map);
distortion_area = area_distortion(v,f,map);

fprintf('Mean(angle distortion) = %.4f\n',mean(abs(distortion_angle)));
fprintf('SD(angle distortion) = %.4f\n',std(abs(distortion_angle)));
fprintf('Mean(area distortion) = %.4f\n',mean(abs(distortion_area)));
fprintf('SD(area distortion) = %.4f\n',std(abs(distortion_area)));

figure;
histogram(distortion_angle,-180:1:180, 'EdgeColor','none'); hold on;
[hist_values, hist_edges] = histcounts(distortion_angle,-180:1:180);
hist_centers = (hist_edges(1:end-1)+hist_edges(2:end))/2;
plot(hist_centers, hist_values)
xlim([-180 180])
title('Angle Distortion');
xlabel('Angle difference (degree)')
ylabel('Number of angles')
set(gca,'FontSize',20);

figure;
histogram(distortion_area,-5:0.1:5, 'EdgeColor','none'); hold on;
[hist_values, hist_edges] = histcounts(distortion_area,-5:0.1:5);
hist_centers = (hist_edges(1:end-1)+hist_edges(2:end))/2;
plot(hist_centers, hist_values)
xlim([-5 5])
title('Area Distortion');
xlabel('log(final area/initial area)')
ylabel('Number of faces')
set(gca,'FontSize',20);

%% Harmonic expansion (Solver)
[thetas, phis] = cart2spheroid(map, foci, zeta, hemispheroid_type);
D_mat = hemispheroidal_harmonic_basis(max_n, thetas, phis, hemispheroid_type);
qm_k = D_mat\new_v;

%% Shape descriptors and the fractal dimension
% Shape descriptors (2-norm) for frequency accumulates at a certain
% frequency degree..
Dl = zeros([3, max_n+1]);
for k_ = 1:3
    for n_ = 1:max_n
        for m_ = -n_:1:n_
            Dl(k_, n_) = Dl(k_, n_) + (real(qm_k(n_^2 + n_ + m_ + 1, k_)))^2 ...
                + (imag(qm_k(n_^2 + n_ + m_ + 1, k_)))^2;
        end
        Dl(k_, n_) = sqrt(Dl(k_, n_));
    end
    Dl(k_, :) = Dl(k_, :)/Dl(k_, 1);
end

figure
subplot(2,3,1)
x_temp = 1:length(Dl(1, :));
loglog(x_temp(2:end), Dl(1, 2:end), 'LineWidth', 2)
title('The shape descriptors')
xlabel('Freqency index (k)')
ylabel('Normalized amplitude (Dx)')
grid on

subplot(2,3,2)
loglog(x_temp(2:end), Dl(2, 2:end), 'LineWidth', 2)
title('The shape descriptors')
xlabel('Freqency index (k)')
ylabel('Normalized amplitude (Dy)')
grid on

subplot(2,3,3)
loglog(x_temp(2:end), Dl(3, 2:end), 'LineWidth', 2)
title('The shape descriptors')
xlabel('Freqency index (k)')
ylabel('Normalized amplitude (Dz)')
grid on

subplot(2,3,4:6)
loglog(x_temp(2:end), sqrt(Dl(1, 2:end).^2 + ...
    Dl(2, 2:end).^2 + Dl(3, 2:end).^2), 'LineWidth', 2)
title('The shape descriptors')
xlabel('Freqency index (k)')
ylabel('Normalized amplitude (Dr)')

%% Harmonic reconstruction
epsillon = pi/75;  % To account for numerical errors.

theta_c = pi /2 - epsillon;
reconstruction_mesh_option = 1;
resolution = 100; % We don't use it with option 1

% Construct reconstruction domain
[rec_v, rec_f] = uniform_spherical_cap_grid(theta_c, edge_length, resolution, reconstruction_mesh_option);
rec_v(:, 3) = rec_v(:, 3) * cc; % This will slightly deform the areas/angles of the tirangles

[rec_thetas, rec_phis] = cart2spheroid(rec_v, foci, zeta, hemispheroid_type);


rec_D_mat = hemispheroidal_harmonic_basis(rec_max_n, rec_thetas, ...
    rec_phis, hemispheroid_type);

rec_v = real(rec_D_mat * qm_k(1:(rec_max_n+1)^2, :));

rec_v_old = rec_v;

% To restore the surface registration
rec_v(:, 3) = normal_sign .* rec_v(:, 3);

rec_v = (R * rec_v')';

rec_v = rec_v + centroid;

paraview_patch(rec_v, rec_f)
fname_out = strcat('./rec_output/rec_n_', num2str(rec_max_n), '_', mapping_type, '_',fname(1:end));
if endsWith(fname_out, '.mat', 'IgnoreCase', true)
    fname_out = strrep(fname_out, '.mat', '.stl');
end
stlWrite(fname_out, rec_f, rec_v);


%% Reconstruct using the original mesh

rec_v_original = real(D_mat * qm_k);
fname_out = strcat('./rec_output/rec_n_orig_', num2str(rec_max_n), '_', mapping_type, '_',fname(1:end));
if endsWith(fname_out, '.mat', 'IgnoreCase', true)
    fname_out = strrep(fname_out, '.mat', '.stl');
end
stlWrite(fname_out, f, rec_v_original);

%% Orthagonality check
degs = linspace(0, (max_n+1)^2, (max_n+1)^2);
[orth_mat1, ~] = meshgrid(degs, degs);
orth_mat2 = orth_mat1;
uniform_D = true;

scatter_vec_x = []; scatter_vec_y = [];

tic
orth_mat1 = real(D_mat' * D_mat);
norms_A = vecnorm(D_mat, 2, 1);
normalization_matrix = norms_A' * norms_A;
orth_mat1 = orth_mat1 ./ normalization_matrix;

orth_mat2 = real(rec_D_mat' * rec_D_mat);
norms_A2 = vecnorm(rec_D_mat, 2, 1);
normalization_matrix2 = norms_A2' * norms_A2;
orth_mat2 = orth_mat2 ./ normalization_matrix2;

toc

viridis_rgb = [
    68, 1, 84;
    72, 35, 116;
    64, 67, 135;
    52, 94, 141;
    41, 120, 142;
    32, 144, 140;
    34, 167, 132;
    68, 190, 112;
    121, 209, 81;
    189, 222, 38;
    253, 231, 37
] / 255;  % Divide by 255 to scale to [0, 1] range

% Interpolate to create a smooth colormap with 256 colors
num_colors = 256;
viridis_colormap = interp1(linspace(0, 1, size(viridis_rgb, 1)), viridis_rgb, linspace(0, 1, num_colors));

figure;
% surf(orth_mat);
imagesc(orth_mat1);
colorbar;
% colormap(parula);
colormap(viridis_colormap);
title('Orthogonality Matrix');
axis equal;
colorbar('FontSize', 16); % Set font size for the colorbar
set(gca, 'FontSize', 16);

hold on
scatter(scatter_vec_x, scatter_vec_y)
hold off

saveas(gcf, 'Mapping_orthogonality.svg');

figure;
% surf(orth_mat);
error_orth = abs(orth_mat1-orth_mat2);
imagesc(error_orth);
colorbar;
% colormap(parula);
colormap(viridis_colormap);
title('Orthogonality Error');
axis equal;
clim([0 1]);
colorbar('FontSize', 16); % Set font size for the colorbar
set(gca, 'FontSize', 16);
saveas(gcf, 'Mapping_orthogonality_Error.svg');

map_face_area = face_area(f, map);
map_face_area = map_face_area ./ sum(map_face_area); % normalized map surface area
area_mean = mean(map_face_area);
area_std = std(map_face_area);

figure;
histo= histogram(map_face_area);
histo.EdgeColor = 'none';

% Customize the histogram appearance (optional)
xlabel('Normalized traingles area', 'FontSize', 16);
ylabel('Frequency', 'FontSize', 16);
title(mapping_type);
set(gca, 'FontSize', 16);
% Hold the current plot to add lines
hold on;

% Add a line for the mean
line([area_mean area_mean], ylim, 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');

% Add lines for the mean ± standard deviation
line([area_mean + area_std area_mean + area_std], ylim, 'Color', 'b', 'LineWidth', 2, 'LineStyle', '--');
line([area_mean - area_std area_mean - area_std], ylim, 'Color', 'b', 'LineWidth', 2, 'LineStyle', '--');

% Add a legend (optional)
legend('Data', 'Mean', 'Mean ± Std Dev', 'Location', 'Best');

% Release the hold on the current plot
hold off;

saveas(gcf, 'Mapping_area_histogram.svg');

rmse = sqrt(mean(error_orth(:).^2));


