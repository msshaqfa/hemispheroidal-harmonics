% Ellipsoidal cap conformal map for simply-connected open surfaces
%
% Usage:
% map = ellipsoidal_cap_conformal_map(v,f,a,b,c)
%
% Input:
% v: nv x 3 vertex coordinates of a simply-connected triangle mesh
% f: nf x 3 triangulations of a simply-connected triangle mesh
% (optional) a,b,c: the radii of the target ellipsoid
%
% Output:
% map: nv x 3 vertex coordinates of the ellipsoidal cap conformal parameterization
%
% Remark:
% - The input surface should be aligned with the xyz axis beforehand,
%   otherwise the result may be affected
%
% Written by Gary Pui-Tung Choi, 2023

addpath('mfile')


%% Example 1: Chinese Lion

load('chinese_lion.mat');

% rotate and align it beforehand
rotx = @(t) [1 0 0; 0 cos(t) -sin(t) ; 0 sin(t) cos(t)] ;
roty = @(t) [cos(t) 0 sin(t) ; 0 1 0 ; -sin(t) 0  cos(t)] ;
rotz = @(t) [cos(t) -sin(t) 0 ; sin(t) cos(t) 0 ; 0 0 1] ;
v = v*rotx(-pi/3.2);
v = bsxfun(@minus,v,mean(v)); % zero-center the data

plot_mesh(v,f,mean_curv); 
view([-15 15]); 
axis on;
xlabel('x')
ylabel('y')
zlabel('z')

% Ellipsoidal cap conformal map
map = ellipsoidal_cap_conformal_map(v,f,1,1,0.9);

% Comparison: Disk conformal map (with Mobius area correction)
% map = disk_conformal_map(v,f);
% map = mobius_area_correction_disk(v,f,map); % with area correction


% plot_mesh(map,f); view([-180 90]);
% can also include the third input if an additional quantity is defined on vertices
plot_mesh(map,f,mean_curv); view([-180 90]); 

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

% rotate and align it beforehand
rotx = @(t) [1 0 0; 0 cos(t) -sin(t) ; 0 sin(t) cos(t)] ;
roty = @(t) [cos(t) 0 sin(t) ; 0 1 0 ; -sin(t) 0  cos(t)] ;
rotz = @(t) [cos(t) -sin(t) 0 ; sin(t) cos(t) 0 ; 0 0 1] ;
v = v*rotx(pi/10)*rotz(pi/2);
v = bsxfun(@minus,v,mean(v)); % zero-center the data

plot_mesh(v,f); view([0 90]);
axis on;
xlabel('x')
ylabel('y')
zlabel('z')

% Ellipsoidal cap conformal map
map = ellipsoidal_cap_conformal_map(v,f,1.1,1,0.9);

% Comparison: Disk conformal map (with Mobius area correction)
% map = disk_conformal_map(v,f);
% map = mobius_area_correction_disk(v,f,map); % with area correction

plot_mesh(map,f); view([-90 90]);

axis on;
xlabel('x')
ylabel('y')
zlabel('z')


%% Example 3: Human Brain

load('human_brain.mat')

% rotate and align it beforehand
rotx = @(t) [1 0 0; 0 cos(t) -sin(t) ; 0 sin(t) cos(t)] ;
roty = @(t) [cos(t) 0 sin(t) ; 0 1 0 ; -sin(t) 0  cos(t)] ;
rotz = @(t) [cos(t) -sin(t) 0 ; sin(t) cos(t) 0 ; 0 0 1] ;
v = v*roty(pi/2);
v = bsxfun(@minus,v,mean(v)); % zero-center the data

% plot_mesh(v,f); view([90 0]);
% can also include the third input if an additional quantity is defined on vertices
plot_mesh(v,f,mean_curv); view([90 0]); 

axis on;
xlabel('x')
ylabel('y')
zlabel('z')

% Ellipsoidal cap conformal map
map = ellipsoidal_cap_conformal_map(v,f,1.5,1,0.75);

% Comparison: Disk conformal map (with Mobius area correction)
% map = disk_conformal_map(v,f);
% map = mobius_area_correction_disk(v,f,map); % with area correction

% plot_mesh(map,f); 
% can also include the third input if an additional quantity is defined on vertices
plot_mesh(map,f,mean_curv); 

axis on;
xlabel('x')
ylabel('y')
zlabel('z')

