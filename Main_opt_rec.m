% Written by Mahmoud Shaqfa and Gary Choi @ 2023
close all; clear; clc;

% addpath('code')
addpath('mfile')
addpath('stlTools')
addpath('distmesh/')
addpath('input_geom/')
addpath('input_geom/others/')
addpath('rec_output/')

%% Load meshes and files

register = true;

 % For loading stl files:
% fname = "3D_face_refined.stl"; % Very wavy reconstruction at only n = 20 ( we need at least n=40 as the SCHA paper)
% fname = "3D_face_refined4.stl"; % Very wavy reconstruction at only n = 20 ( we need at least n=40 as the SCHA paper)
% fname = "Matterhorn_new_mode.stl"; % Very wavy reconstruction at n = 20 (we need n=60 as DHA paper)

fname = "half_stone.stl"; % The reconstruction is really good here (n = 30)!
% alpha = 0.9306, beta = 0.0693, gamma = 1.5882e-04 ( n = 15)
% alpha = 0.4795, beta = 0.4795, gamma = 0.0410 ( n = 5)
% 0.5647    0.4313    0.0041 (n = 20)


% fname = "Ziad_stone_patch2.stl"; % fails for rough surfaces in general
% fname = "3D_face_refined4_filled_v3.stl";
% fname = "half_human_femur.stl";
% fname = "human_face.stl";

[v, f, ~, ~] = stlRead(fname); % Uncomment only when loading STL files.



% For loading *.mat files
% fname = 'chinese_lion.mat'; % The reconstruction kind of acceptable (tested to n = 30)!
% fname = 'human_brain.mat'; % Completely fails in reconstruction
% fname = 'human_face.mat'; % The nose reconstruction will be very complicated (fails)

% load(fname);  % Uncomment only when loading MAT files.



% Parameter setting
max_n = 30; % Max Analysis degree
rec_max_n = max_n; % Max reconstruction degree

opt_rec_max_n = 10; % For the optimization problem

edge_length = 0.0350; % Reconstruction domain resolution

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

% err
%% Fit hemispheroidal shape over the mesh
% This part needs more work (Discuss with Gary the optimization)

% Bounding box after alignment
BB_x = max(new_v(:, 1)) - min(new_v(:, 1));
BB_y = max(new_v(:, 2)) - min(new_v(:, 2));
BB_z = max(new_v(:, 3)) - min(new_v(:, 3));

BB_xy = mean([BB_x, BB_y]);
BB_z_normalized = BB_z / BB_xy;

if BB_z_normalized >= 1
    hemispheroid_type = "prolate";
else
    hemispheroid_type = "oblate";
end

%% Compute the spheroidal coordinates parameters (with a unit disk base)
foci = sqrt(abs(1 - BB_z_normalized^2)); % Focal distance
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
%% Optimize ellipsoidal mapping (Optimal reconstruction via balance mapping)


% Precompute mappings
map_T = hemispheroidal_tutte_map(new_v, f, cc);
map_C = hemispheroidal_conformal_map(new_v, f, cc);        
map_A = hemispheroidal_area_preserving_map(new_v, f, cc);

hemispheroid_data = [];
hemispheroid_data.foci = foci;
hemispheroid_data.zeta = zeta;
hemispheroid_data.hemispheroid_type = hemispheroid_type;
hemispheroid_data.c = cc; % This is redundant.


% Optimize balanced mapping
disp("The optimization stage")
tic
[map,alpha,beta,gamma,qm_k_i] = opt_map4rec(v, f, map_T, map_C, map_A, opt_rec_max_n, hemispheroid_data);
toc

disp('Optimal alpha, beta, gamma:')
disp([alpha, beta, gamma])

%% Harmonic expansion (Solver)
[thetas, phis] = cart2spheroid(map, foci, zeta, hemispheroid_type);
D_mat = hemispheroidal_harmonic_basis(max_n, thetas, phis, hemispheroid_type, false);
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
epsillon = pi/30;  % To account for numerical errors.
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

% To restore the surface registration
rec_v(:, 3) = normal_sign .* rec_v(:, 3);

rec_v = (R * rec_v')';
% rec_v = (rec_v * R');

rec_v = rec_v + centroid;
paraview_patch(rec_v, rec_f)
fname_out = strcat('./rec_output/opt_rec_', '_',fname(1:end));
if endsWith(fname_out, '.mat', 'IgnoreCase', true)
    fname_out = strrep(fname_out, '.mat', '.stl');
end
stlWrite(fname_out, rec_f, rec_v);

