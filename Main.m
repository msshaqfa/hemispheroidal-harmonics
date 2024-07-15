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
% fname = "3D_face_refined.stl"; % Very wavy reconstruction at only n = 20 ( we need at least n=40 as the SCHA paper)
% fname = "3D_face_refined4.stl"; % Very wavy reconstruction at only n = 20 ( we need at least n=40 as the SCHA paper)
fname = "Matterhorn_new_mode.stl"; % Very wavy reconstruction at n = 20 (we need n=60 as DHA paper)
% fname = "half_stone.stl"; % The reconstruction is really good here (n = 30)!
% fname = "Ziad_stone_patch2.stl"; % fails for rough surfaces in general
% fname = "3D_face_refined4_filled_v3.stl";
% fname = "half_human_femur.stl";
% fname = "human_face.stl"; % The input mesh is very sparse and reconstruction results are bad
% fname = "bumpy_ref_AR_3_open.stl";
% fname = "3D_face_simplified.stl";
% fname = "sophie.stl";
% fname = "bunny_open.stl";

[v, f, ~, ~] = stlRead(fname); % Uncomment only when loading STL files.

% For loading *.mat files
% fname = 'chinese_lion.mat'; % The reconstruction kind of acceptable (tested to n = 30)!
% fname = 'human_brain.mat'; % Completely fails in reconstruction
% fname = 'human_face.mat'; % The nose reconstruction will be very complicated (fails)

% Teeth benchmarks
% fname = "a10_sas_aligned.mat";
% fname = "D09_sas_aligned.mat";
% fname = "S09_sas_aligned.mat"; % very good
% fname = "T12_sas_aligned.mat"; % very good
% fname = "u16_sas_aligned.mat"; % Perfect with AP
% fname = "V09_sas_aligned.mat"; % good
% fname = "w02_sas_aligned.mat"; % good
% fname = "x02_sas_aligned.mat"; % good
% fname = "t09_sas_aligned.mat"; % good
% fname = "Q12_sas_aligned.mat"; % good
% fname = "P32_sas_aligned.mat";
% fname = "k18_sas_aligned.mat";
% fname = "J12_sas_aligned.mat";
% fname = "i23_sas_aligned.mat"; % all good until here
% fname = "H19_sas_aligned.mat"; % AP is magic
% fname = "B03_sas_aligned.mat"; % good maybe alignement need to be fiexd


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
% alpha = 0.4071;
alpha = 0.35;
beta = 0.0;
gamma = 0.65;

% optimal for bunny:
% beta = 0.3818;
% gamma = 0.6182;

% optimal for 3D_face
% alpha = 0.35;
% gamma = 0.65;


% optimal for 3D_face (RMSE=0.00149)
% alpha = 0.40;
% gamma = 0.60;

% optimal for 3D_face (RMSE=0.00156)
% alpha = 0.45;
% gamma = 0.65;

% Optimal mapping: 0.3282 0 0.6718 (RMS : 0.001553)
% Optimal mapping: 0.3440 0 0.6560 (RMS : 0.001544)
% Optimal mapping: 0.5 0.5 0.0 (RMS : )

% For Matterhorn 
% optimal mapping: 0.3377    0.3377    0.3245 ()
% optimal mapping: 0.3333    0.0       0.6667 ()

% Parameter setting
max_n = 10; % Max Analysis degree
rec_max_n = max_n; % Max reconstruction degree
% edge_length = 0.025; % Reconstruction domain resolution
% edge_length = 0.0170; % Reconstruction domain resolution (bunny)
edge_length = 0.0250;% Reconstruction domain resolution (face)
% edge_length = 0.010;%250; % Reconstruction domain resolution (Brain)

% plot_mesh(v,f,mean_curv); 
plot_mesh(v,f);
view([-15 15]); 
axis on;
xlabel('x')
ylabel('y')
zlabel('z')




%% Transfer data to cononical coordinates (coordinates shift and rotation)
% Surface registration step:
[new_v, centroid, R, normal_sign] = register_surface(v, f, true);
% To restore the rotation of the coordinates use: (R * new_v)'

% if normal_sign < 0
%     new_v(:, 3) = new_v(:, 3) * normal_sign;
% end

if max(new_v(:, 3)) ~= max(abs(new_v(:, 3)))
    new_v(:, 3) = new_v(:, 3) * -1;
end

% plot_mesh(new_v, f, mean_curv);
plot_mesh(new_v, f);

view([-15 15]); 
axis on;
xlabel('x')
ylabel('y')
zlabel('z')


% 
% [~, outer] = meshboundaries(f);
% 
% fig = figure;
% subplot(1,2,1);
% % patch('Faces',f,'Vertices',v,'FaceColor',[200 191 231]/255,'EdgeColor','none');
% h1 = patch('Faces',f,'Vertices',v,'FaceColor',[200 191 231]/255,'EdgeColor','black');
% 
% % hold on;
% % bds_v = v(outer, :);
% % plot3(bds_v(:,1), bds_v(:,2), bds_v(:,3), 'b', 'LineWidth', 4);
% % hold off;
% % axis equal tight
% axis equal tight off 
% 
% % view([-100 10]) % for brain
% view([-21.1000, 24.1389]) % for the face
% xlabel("x")
% ylabel("y")
% zlabel("z")
% box on;
% camlight;
% % Increase line width of axes
% set(gca, 'LineWidth', 1);
% 
% % Change font size of ticks and labels
% set(gca, 'FontSize', 14); % Set font size for tick labels
% 
% 
% subplot(1,2,2);
% % patch('Faces',f,'Vertices',new_v,'FaceColor',[200 191 231]/255,'EdgeColor','none');
% patch('Faces',f,'Vertices',new_v,'FaceColor',[200 191 231]/255,'EdgeColor','black');
% % hold on;
% % bds_v = new_v(outer, :);
% % plot3(bds_v(:,1), bds_v(:,2), bds_v(:,3), 'r', 'LineWidth', 4);
% % hold off;
% % axis equal tight
% axis equal tight off 
% % view([-100 10]) % for brain
% % [current_azimuth, current_elevation] = view;
% view([-21.1000, 24.1389])
% xlabel("x")
% ylabel("y")
% zlabel("z")
% box on;
% camlight;
% 
% % Increase line width of axes
% set(gca, 'LineWidth', 2);
% set(fig, 'Color', 'none');
% 
% % Change font size of ticks and labels
% set(gca, 'FontSize', 14); % Set font size for tick labels
% print('-dsvg', 'surfaces_vis.svg', '-r300', '-S500,400', '-noui', '-painters', '-opengl');

% err
%% Fit hemispheroidal shape over the mesh
% This part needs more work (Discuss with Gary the optimization)

% Bounding box after alignment
BB_d = max(new_v) - min(new_v);
BB_x = BB_d(1); BB_y = BB_d(2); BB_z = BB_d(3);
% BB_x = max(new_v(:, 1)) - min(new_v(:, 1));
% BB_y = max(new_v(:, 2)) - min(new_v(:, 2));
% BB_z = max(new_v(:, 3)) - min(new_v(:, 3)); % height of the BB

% BB_z = mean(new_v(:, 3)); % mean height of the BB

BB_xy = mean([BB_x, BB_y]);
BB_z_normalized = abs(BB_z) / BB_xy;

% For tests to be deleted:
BB_z_normalized = 0.9999;

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
    % title("Fitted oblate hemispheroid")
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
epsillon = pi/75;  % To account for numerical errors. (Matterhorn)

% epsillon = pi/160;  % To account for numerical errors. (3D Face via AP & Tutte)
epsillon = pi/100;  % To account for numerical errors. (3D Face via conformal)

% epsillon = pi/20;  

theta_c = pi /2 - epsillon;
reconstruction_mesh_option = 1;
resolution = 100; % We don't use it with option 1

% rec_max_n = 2; % re-write for testing
% rec_max_n = 20;

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
% rec_v = (rec_v * R');

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

%% Render figures

% err
% fname_in = strcat('./rec_output/rec_n_', num2str(100), '_', mapping_type, '_',fname(1:end));
% fname_in = 'rec_n_75_balanced_3D_face_refined.stl';
% fname_in = 'rec_n_75_tutte_3D_face_refined.stl';
% fname_in = 'rec_n_40_conformal_3D_face_refined.stl';
% fname_in = 'rec_n_40_area_preserving_3D_face_refined.stl';
fname_in = 'rec_n_50_area_preserving_Matterhorn_new_mode.stl';
fname_in = 'rec_n_50_area_preserving_Matterhorn_new_mode_0.6250.stl';

% fname_in = 'rec_n_75_tutte_3D_face_refined.stl';
% fname_in = 'rec_n_50_balanced_Matterhorn_new_mode.stl';

% fname_in = 'rec_n_100_area_preserving_human_brain.stl';
[rec_v_old, rec_f, ~, ~] = stlRead(fname_in);
rec_v_old = (rec_v_old * R);

rec_v_old(:, 3) = -rec_v_old(:, 3); % For Matterhorn example.
fig = figure;
% subplot(1,2,1);

h1 = patch('Faces',rec_f,'Vertices',rec_v_old,'FaceColor',[200 191 231]/255,'EdgeColor','black', 'LineWidth', 0.5);
% h1 = patch('Faces',rec_f,'Vertices',rec_v_old,'FaceColor',[200 191 231]/255,'EdgeColor','none');

% hold on;
% bds_v = v(outer, :);
% plot3(bds_v(:,1), bds_v(:,2), bds_v(:,3), 'b', 'LineWidth', 4);
% hold off;
% axis equal tight
axis equal tight off 

% view([-10 90]) % for brain
% view([-90 -75]) % for 3D face
view([-525 16]) % for Matterhorn
% view([220 15]) % for bunny

% view([155 35]) % for bunny 2
% view([115 30]) % for bunny 3

xlabel("x")
ylabel("y")
zlabel("z")
box on;
camlight;
% Increase line width of axes
% set(gca, 'LineWidth', 0.5);

% Change font size of ticks and labels
set(gca, 'FontSize', 14); % Set font size for tick labels





%% Orthagonality
% err
degs = linspace(0, (max_n+1)^2, (max_n+1)^2);
[orth_mat1, ~] = meshgrid(degs, degs);
orth_mat2 = orth_mat1;
uniform_D = true;

scatter_vec_x = []; scatter_vec_y = [];
% tic
% % if ~uniform_D
%     for n1 = 0 : max_n
%         for n2 = 0 : max_n
%             for m1 = -n1:n1
%                 for m2 = -n2:n2
%                     id1 = n1 * n1 + n1 + m1 + 1;
%                     id2 = n2 * n2 + n2 + m2 +1;
%                     norm1 = norm(D_mat(:, id1));
%                     norm2 = norm(D_mat(:, id2));
%                     orth_mat1(id1, id2) = real(D_mat(:, id1)' * D_mat(:, id2) / (norm1 * norm2));
%                     % orth_mat(id1, id2) = real(D_mat(:, id1)' * D_mat(:, id2));
%                     if (m1==0 && m2==0)
%                         idx = length(scatter_vec_y);
%                         scatter_vec_x(idx+1) = id1;
%                         scatter_vec_y(idx+1) = id2;
%                     end
%                 end
%             end
%         end
%     end
%     toc
% else
    % for n1 = 0 : max_n
    %     for n2 = 0 : max_n
    %         for m1 = -n1:n1
    %             for m2 = -n2:n2
    %                 id1 = n1 * n1 + n1 + m1 + 1;
    %                 id2 = n2 * n2 + n2 + m2 +1;
    %                 norm1 = norm(rec_D_mat(:, id1));
    %                 norm2 = norm(rec_D_mat(:, id2));
    %                 orth_mat2(id1, id2) = real(rec_D_mat(:, id1)' * rec_D_mat(:, id2) / (norm1 * norm2));
    %                 % orth_mat(id1, id2) = real(D_mat(:, id1)' * D_mat(:, id2));         
    %             end
    %         end
    %     end
    % end
    % 
% end

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
% 
% hold on;
% rectangles = [];
% for i = 0:max_n
%     val = (i+1)^2;
%     pos = ( i )^2;
%     plot([val, 0], [0, val], 'LineWidth', 5);
% end
% 
% hold off;

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

% grid on;

rmse = sqrt(mean(error_orth(:).^2))


xx = [0.1827
0.3653
0.4
0.5
0.625
0.75
0.85
0.9999];

yy = [0.000662
0.000622
0.000578
0.000572
0.000561
0.00054
0.00054
0.00057];

yy2 = [0.0332
0.0218
0.0204
0.0177
0.0171
0.0186
0.0204
0.0231];

figure;

% Plot the first data set with respect to the left y-axis
yyaxis left;
semilogy(xx, yy, '-b', 'DisplayName', 'A-RMSE');
ylabel('A-RMSE');
ylim([min(yy), max(yy)]);
set(gca, 'FontSize', 16);
% grid on;

% figure;
% Plot the second data set with respect to the right y-axis
yyaxis right;
semilogy(xx, yy2, '-r', 'DisplayName', 'RMSE');
ylabel('RMSE');
ylim([min(yy2), max(yy2)]);

% Add labels and title
xlabel('Hemispheroidal depth c');
% set(gca, 'FontSize', 16);
% grid on;
% Add a legend
legend show;



