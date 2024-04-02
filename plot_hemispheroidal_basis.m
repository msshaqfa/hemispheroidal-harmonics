% Mahmoud Shaqfa

clear; close all; clc;


n = 8;
m = 4;
grid = 300;
% grid = 70;

a = 5;
zeta = 0.4;
e = 0.5;

% Create the grid
theta= linspace(0, pi/2, grid);
phi = linspace(0, 2.0 * pi - 0.0001, grid);
[phi, theta] = meshgrid(phi, theta);


%% Colormaps (paraview)
% Define a color map similar to Paraview's shades
red_color = zeros(1, 255); green_color = zeros(1, 255); blue_color = zeros(1, 255);
red_color(1:127) = linspace(42,220,127); red_color(128:255) = linspace(220,174,127+1);
green_color(1:127) = linspace(63,220,127); green_color(128:255) = linspace(220,0,127+1);
blue_color(1:127) = linspace(181,220,127); blue_color(128:255) = linspace(220,22,127+1);
ParaviewMap = [red_color', green_color', blue_color']./255;



%% Oblate hemispheroid harmonics

% Draw surface
x = e .* cosh(zeta) .* cos(theta) .* cos(phi);
y = e .* cosh(zeta) .* cos(theta) .* sin(phi);
z = e .* sinh(zeta) .* sin(theta);


% Harmonic function
N_mn = sqrt((2*n +1) * factorial(n-m) / (4*pi * factorial(n+m)));
P_mn = legendre(n, 2.* sin(theta) - 1);
P_mn = N_mn .* real(squeeze(P_mn(m+1, :, :)) .* exp(m .* 1i .* phi));


% Plot the surface
figure(Renderer="painters");
h1 = surf(x,y,z, P_mn);
set(h1,'edgecolor','none')


% hold on;
% surf(x,y,-z, P_mn)

light
lighting phong
axis tight equal off
view(40,30)
camzoom(1.3)
colormap(ParaviewMap);
saveas(gcf, strcat('oblate_n_', num2str(n), '_m_', num2str(m), '.svg'));


%% Prolate harmonics

theta2 = linspace(0, .5 * pi - 0.0001, grid);
phi2 = linspace(0, 2.0 * pi - 0.0001, grid);
[phi2, theta2] = meshgrid(phi2, theta2);

% Draw surface
x2 = e .* sinh(zeta) .* sin(theta2) .* cos(phi2);
y2 = e .* sinh(zeta) .* sin(theta2) .* sin(phi2);
z2 = e .* cosh(zeta) .* cos(theta2);


% Harmonic function
N_mn2 = sqrt((2*n +1) * factorial(n-m) / (4*pi * factorial(n+m)));
P_mn2 = legendre(n, 1 - 2.* cos(theta2));
P_mn2 = N_mn2 .* real(squeeze(P_mn2(m+1, :, :)) .* exp(m .* 1i .* phi2));


% Plot the surface
figure(Renderer="painters");
h2 = surf(x2,y2,z2, P_mn2);
set(h2,'edgecolor','none')

light
lighting phong
axis tight equal off
view(40,30)
camzoom(1.3)
colormap(ParaviewMap);
saveas(gcf, strcat('prolate_n_', num2str(n), '_m_', num2str(m), '.svg'));




%% Verification for Wim (simple stretched/compressed spherical harmonics)

% For oblate
scale_xy_o = e .* cosh(zeta);
scale_z_o = e .* sinh(zeta);

scale_xy_p = e .* sinh(zeta);
scale_z_p = e .* cosh(zeta);


x_s = cos(theta) .* cos(phi);
y_s = cos(theta) .* sin(phi);
z_s = sin(theta);




subplot(3, 2, 1)
h1 = surf(x,y,z, P_mn);
set(h1,'edgecolor','none')

light
lighting phong
axis tight equal off
view(40,30)
camzoom(1.3)
colormap(ParaviewMap);
title("Oblate")

subplot(3, 2, 2)
h2 = surf(x2,y2,z2, P_mn2);
set(h2,'edgecolor','none')

light
lighting phong
axis tight equal off
view(40,30)
camzoom(1.3)
colormap(ParaviewMap);
title("Prolate")



N_mn = sqrt((2*n +1) * factorial(n-m) / (4*pi * factorial(n+m)));
P_mn_hemi = legendre(n, 1 - 2.* cos(theta));
P_mn_hemi = N_mn .* real(squeeze(P_mn_hemi(m+1, :, :)) .* exp(m .* 1i .* phi));


subplot(3, 2, 3:4)
h3 = surf(x_s,y_s,z_s, P_mn_hemi);
set(h3,'edgecolor','none')

light
lighting phong
axis tight equal off
view(40,30)
camzoom(1.3)
colormap(ParaviewMap);
title("Hemisphere")




subplot(3, 2, 5)
h4 = surf(scale_xy_o .* x_s, scale_xy_o .* y_s, scale_z_o .* z_s, P_mn_hemi);
set(h4,'edgecolor','none')

light
lighting phong
axis tight equal off
view(40,30)
camzoom(1.3)
colormap(ParaviewMap);
title("Scaled Oblate Hemisphere")





subplot(3, 2, 6)
h5 = surf(scale_xy_p .* x_s, scale_xy_p .* y_s, scale_z_p .* z_s, P_mn_hemi);
set(h5,'edgecolor','none')

light
lighting phong
axis tight equal off
view(40,30)
camzoom(1.3)
colormap(ParaviewMap);
title("Scaled Prolate Hemisphere")





