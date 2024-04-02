% Written by Mahmoud Shaqfa and Gary Choi @ 2023
close all; clear; clc;

% addpath('code')
addpath('mfile')
addpath('stlTools')
addpath('distmesh/')
addpath('input_geom/')
addpath('input_geom/others/')
addpath('rec_output/')
addpath('input_geom/Gary/') % Gary Tests

%% Load meshes and files

register = false;

 % For loading stl files:

% 3D face example
% fname = "3D_face_refined4_filled_clean.off";
% fname_map = "3D_face_refined4_filled_tutte_map.off";
% fname_map = "3D_face_refined4_filled_conformal_map.off";
% fname_map = "3D_face_refined4_filled_area_preserving_map.off";


% Half-stone example
% fname = "half_stone_clean.off";
% fname_map = "half_stone_tutte_map.off";
% fname_map = "half_stone_conformal_map.off";
% fname_map = "half_stone_area_preserving_map.off";

% Matterhorn example
fname = "Matterhorn_new_mode_clean.off";
% fname_map = "Matterhorn_new_mode_tutte_map.off";
% fname_map = "Matterhorn_new_mode_conformal_map.off";
fname_map = "Matterhorn_new_mode_area_preserving_map.off";



% [v, f, ~, ~] = stlRead(fname); % Uncomment only when loading STL files.
% [map, ~, ~, ~] = stlRead(fname_map); % Uncomment only when loading STL files.

[v, f] = read_off(fname); % Uncomment only when loading STL files.
v = v';
f = f';
[map, ~] = read_off(fname_map); % Uncomment only when loading STL files.
map = map';

% For loading *.mat files
% fname = 'chinese_lion.mat'; % The reconstruction kind of acceptable (tested to n = 30)!
% fname = 'human_brain.mat'; % Completely fails in reconstruction
% fname = 'human_face.mat'; % The nose reconstruction will be very complicated (fails)

% Teeth benchmarks
% fname = "a10_sas_aligned.mat";
% fname = "D09_sas_aligned.mat";
% fname = "S09_sas_aligned.mat"; % very good
% fname = "T12_sas_aligned.mat"; % very good
% fname = "u16_sas_aligned.mat"; % kinda fails
% fname = "V09_sas_aligned.mat"; % good
% fname = "w02_sas_aligned.mat"; % good
% fname = "x02_sas_aligned.mat"; % good
% fname = "t09_sas_aligned.mat"; % good
% fname = "Q12_sas_aligned.mat"; % good
% fname = "P32_sas_aligned.mat";
% fname = "k18_sas_aligned.mat";
% fname = "J12_sas_aligned.mat";
% fname = "i23_sas_aligned.mat"; % all good until here
% fname = "H19_sas_aligned.mat"; % fails
% fname = "B03_sas_aligned.mat"; % good maybe alignement need to be fiexd


% load(fname);  % Uncomment only when loading MAT files.



% Parameter setting
max_n = 40; % Max Analysis degree
rec_max_n = max_n; % Max reconstruction degree
edge_length = 0.0350; % Reconstruction domain resolution
edge_length = 0.0250; % Reconstruction domain resolution

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

if normal_sign < 0
    new_v(:, 3) = new_v(:, 3) * normal_sign;
end

% plot_mesh(new_v, f, mean_curv);
plot_mesh(new_v, f);

view([-15 15]); 
axis on;
xlabel('x')
ylabel('y')
zlabel('z')


%% Fit hemispheroidal shape over the mesh
% This part needs more work (Discuss with Gary the optimization)

% Bounding box after alignment
aa = (max(map(:, 1)) - min(map(:, 1))) / 2;
bb = (max(map(:, 2)) - min(map(:, 2))) / 2;
cc = max(map(:, 3)) - min(map(:, 3));
 
if aa > cc
    hemispheroid_type = "oblate";
else
    hemispheroid_type = "prolate";
end

%% Compute the spheroidal coordinates parameters (with a unit disk base)
% foci = sqrt(abs(1 - BB_z_normalized^2)); % Focal distance
foci = sqrt(abs(aa^2 - cc^2)); % Focal distance
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
    title("Fitted oblate hemispheroid")
    axis equal
end
disp(strcat("The surface is ", hemispheroid_type, " with a = ",...
    num2str(aa), " and c = ", num2str(cc)))
%% Ellipsoidal cap conformal map


% fname_map_out = strcat('./rec_output/map_', fname(1:end));
% stlWrite(fname_map_out, f, map);

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
histogram(distortion_angle,-180:1:180); hold on;
[hist_values, hist_edges] = histcounts(distortion_angle,-180:1:180);
hist_centers = (hist_edges(1:end-1)+hist_edges(2:end))/2;
plot(hist_centers, hist_values)
xlim([-180 180])
title('Angle Distortion');
xlabel('Angle difference (degree)')
ylabel('Number of angles')
set(gca,'FontSize',12);

figure;
histogram(distortion_area,-5:0.1:5); hold on;
[hist_values, hist_edges] = histcounts(distortion_area,-5:0.1:5);
hist_centers = (hist_edges(1:end-1)+hist_edges(2:end))/2;
plot(hist_centers, hist_values)
xlim([-5 5])
title('Area Distortion');
xlabel('log(final area/initial area)')
ylabel('Number of faces')
set(gca,'FontSize',12);

% err
%% Harmonic expansion (Solver)
[thetas, phis] = cart2spheroid(map, foci, zeta, hemispheroid_type);
D_mat = hemispheroidal_harmonic_basis(max_n, thetas, phis, hemispheroid_type);
% qm_k = D_mat\new_v;
qm_k = D_mat\v;

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
epsillon = pi/70;  % To account for numerical errors.
theta_c = pi /2 - epsillon;
reconstruction_mesh_option = 1;
resolution = 100; % We don't use it with option 1

% rec_max_n = 2; % re-write for testing

% Construct reconstruction domain
[rec_v, rec_f] = uniform_spherical_cap_grid(theta_c, edge_length, resolution, reconstruction_mesh_option);
rec_v(:, 3) = rec_v(:, 3) * cc; % This will slightly deform the areas/angles of the tirangles

[rec_thetas, rec_phis] = cart2spheroid(rec_v, foci, zeta, hemispheroid_type);

rec_D_mat = hemispheroidal_harmonic_basis(rec_max_n, rec_thetas, ...
    rec_phis, hemispheroid_type);

rec_v = real(rec_D_mat * qm_k(1:(rec_max_n+1)^2, :));
% To restore the rotation of the coordinates use: (R * new_v)'

paraview_patch(rec_v, rec_f)
fname_out = strcat('./rec_output/rec_', fname_map(1:end));

write_off(fname_out, rec_v', rec_f');





