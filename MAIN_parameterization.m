% Test program for the proposed hemispheroidal parameterization methods:
% - Tutte: hemispheroidal_tutte_map.m
% - Conformal: hemispheroidal_conformal_map.m
% - Area-preserving: hemispheroidal_area_preserving_map.m
% - Balanced: hemispheroidal_balanced_map.m
%
% Written by Gary Pui-Tung Choi, 2024

addpath('mfile')
addpath(genpath('input_geom'));
addpath(genpath('io'));

%% Load a mesh

load('example_face.mat')
% load('example_sophie.mat')
% load('example_julius.mat')
% load('example_chinese_lion.mat')
% load('example_matterhorn.mat')
% load('example_brain.mat')
% load('example_bunny.mat')


% clean the mesh and remove all valence 2 vertices (vertices connected to only 1 triangle) at the boundary
[v,f] = clean_mesh(v,f,1); 


plot_mesh(v,f); 

%% Determine a suitable hemispheroidal shape

% Transfer data to cononical coordinates (coordinates shift and rotation)
% Surface registration step:
[new_v, centroid, R, normal_sign] = register_surface(v, f, false);
% To restore the rotation of the coordinates use: (R * new_v)'

if normal_sign < 0
    new_v(:, 3) = new_v(:, 3) * normal_sign;
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

