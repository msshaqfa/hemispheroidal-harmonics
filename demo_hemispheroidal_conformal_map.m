% Hemispheroidal conformal map for simply-connected open surfaces
%
% Usage:
% map = hemispheroidal_conformal_map(v,f,c)
%
% Input:
% v: nv x 3 vertex coordinates of a simply-connected triangle mesh
% f: nf x 3 triangulations of a simply-connected triangle mesh
% c: the height of the target hemispheroid
%
% Output:
% map: nv x 3 vertex coordinates of the hemispheroidal conformal parameterization
%
% Remark:
% - The input surface should be aligned with the xyz axis beforehand,
%   otherwise the result may be affected
%
% Written by Gary Pui-Tung Choi, 2023

addpath('mfile')


%% Example 1: Chinese Lion

load('chinese_lion.mat');

plot_mesh(v,f,mean_curv); 
view([-15 15]); 
axis on;
xlabel('x')
ylabel('y')
zlabel('z')

% Hemispheroidal conformal map
map = hemispheroidal_conformal_map(v,f,1.5);

% plot_mesh(map,f); view([-180 90]);
% can also include the third input if an additional quantity is defined on vertices
plot_mesh(map,f,mean_curv); view([-180 90]); 

axis on;
xlabel('x')
ylabel('y')
zlabel('z')
view([150 20]);

%% Evaluate the angle and area distortion 
distortion_angle = angle_distortion(v,f,map);
distortion_area = area_distortion(v,f,map);

fprintf('Mean(angle distortion) = %.4f\n',mean(abs(distortion_angle)));
fprintf('SD(angle distortion) = %.4f\n',std(abs(distortion_angle)));
fprintf('Mean(area distortion) = %.4f\n',mean(abs(distortion_area)));
fprintf('SD(area distortion) = %.4f\n',std(abs(distortion_area)));

figure;
histogram(distortion_angle,-180:1:180);
xlim([-180 180])
title('Angle Distortion');
xlabel('Angle difference (degree)')
ylabel('Number of angles')
set(gca,'FontSize',12);

figure;
histogram(distortion_area,-5:0.1:5);
xlim([-5 5])
title('Area Distortion');
xlabel('log(final area/initial area)')
ylabel('Number of faces')
set(gca,'FontSize',12);


%% Example 2: Face

load('human_face.mat');

plot_mesh(v,f); view([0 90]);
axis on;
xlabel('x')
ylabel('y')
zlabel('z')

% Hemispheroidal conformal map
map = hemispheroidal_conformal_map(v,f,0.8);

plot_mesh(map,f); view([-90 90]);

axis on;
xlabel('x')
ylabel('y')
zlabel('z')
view([150 20]);

%% Example 3: Human Brain

load('human_brain.mat')

% plot_mesh(v,f); view([90 0]);
% can also include the third input if an additional quantity is defined on vertices
plot_mesh(v,f,mean_curv); view([90 0]); 

axis on;
xlabel('x')
ylabel('y')
zlabel('z')

% Hemispheroidal conformal map
map = hemispheroidal_conformal_map(v,f,2);

% plot_mesh(map,f); 
% can also include the third input if an additional quantity is defined on vertices
plot_mesh(map,f,mean_curv); 

axis on;
xlabel('x')
ylabel('y')
zlabel('z')
view([150 20]);

