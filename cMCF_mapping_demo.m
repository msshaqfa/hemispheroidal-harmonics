% Written by Mahmoud Shaqfa and Gary Choi @ 2023
close all; clear; clc;

% addpath('code')
addpath('mfile')
addpath('stlTools')
addpath('distmesh/')
addpath('input_geom/')
addpath('rec_output/')




%% Load meshes and files

register = true;

 % For loading stl files:
% fname = "3D_face_refined.stl"; % Very wavy reconstruction at only n = 20 ( we need at least n=40 as the SCHA paper)
fname = "Matterhorn_new_mode.stl"; % Very wavy reconstruction at n = 20 (we need n=60 as DHA paper)
% fname = "half_stone.stl"; % The reconstruction is really good here (n = 30)!

[v, f, ~, ~] = stlRead(fname); % Uncomment only when loading STL files.


% For loading *.mat files
% fname = 'chinese_lion.mat'; % The reconstruction kind of acceptable (tested to n = 30)!
% fname = 'human_brain.mat'; % Completely fails in reconstruction
% fname = 'human_face.mat'; % The nose reconstruction will be very complicated

% load(fname);  % Uncomment only when loading MAT files.



% Parameter setting
max_n = 20; % Max Analysis degree
rec_max_n = max_n; % Max reconstruction degree
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
title("Original mesh before cMCF")
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


%% Conformalized Mena Curvature Flow (cMCF) approach

iterations = 5;
time_step = 0.0005;

M = lumped_mass_matrix(new_v, f);   % Mass matrix (Diagonalized)
Minv = 1 \ M;
L = cotangent_laplacian(new_v, f)/4.0;  % Cotangent Laplacian operator (we can try other defs.)

% spy(M)
% title("Diagonalized mass matrix (M)")
% spy(L)
% title("Cot-Laplacian matrix (L)")

U = new_v;
original_area = sum(face_area(f, U));
original_scale = sqrt(original_area); % Surface scale

original_centroid = mean(U); % Centroid of boundaries
new_v = new_v - original_centroid;


for i = 1:iterations
    M = lumped_mass_matrix(U, f);   % Mass matrix (Diagonalized)
    % Solve (M-delta*L) U = M*U
    S = (M  - time_step .* L);
    U = S \ (M * U);
    
    centroid = mean(U);
    
    % Normalize for the centroid
    U = U - centroid;


    % Normalize the surface area
    area = sum(face_area(f, U));
    disp(sqrt(area))
    
    U = U ./ sqrt(area);

end

U_ = real(U);

% Rescale the mesh
U_map = U_ .* original_scale + centroid;


% plot_mesh(new_v, f, mean_curv);
plot_mesh(U_map, f);
title("Refined mesh after cMCF")
view([-15 15]); 
axis on;
xlabel('x')
ylabel('y')
zlabel('z')



% Evaluate the angle and area distortion 
distortion_angle = real(angle_distortion(new_v, f, U_map));
distortion_area = real(area_distortion(new_v, f, U_map));

figure;
subplot(1, 2, 1)
histogram(distortion_angle,-180:1:180); hold on;
[hist_values, hist_edges] = histcounts(distortion_angle,-180:1:180);
hist_centers = (hist_edges(1:end-1)+hist_edges(2:end))/2;
plot(hist_centers, hist_values)
xlim([-180 180])
title('Angle Distortion');
xlabel('Angle difference (degree)')
ylabel('Number of angles')
set(gca,'FontSize',12);


subplot(1, 2, 2)
histogram(distortion_area,-5:0.1:5); hold on;
[hist_values, hist_edges] = histcounts(distortion_area,-5:0.1:5);
hist_centers = (hist_edges(1:end-1)+hist_edges(2:end))/2;
plot(hist_centers, hist_values)
xlim([-5 5])
title('Area Distortion');
xlabel('log(final area/initial area)')
ylabel('Number of faces')
set(gca,'FontSize',12);


