function map = ellipsoidal_cap_conformal_map(v,f,a,b,c)
% A fast method for computing ellipsoidal cap conformal map of an open surface
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

if nargin < 5
    % Optimally rotate the input surface
    v_ori = v;
    A = bsxfun(@minus,v_ori,mean(v_ori)); % center the data
    [~,~,V] = svd(A,0);
    v = (A*V)*[0,0,1;0,1,0;-1,0,0]; % first align with x-axis then z-axis

    % Use a rectangular bounding box to define a,b,c
    a = (max(v(:,1))-min(v(:,1)))/2;
    b = (max(v(:,2))-min(v(:,2)))/2;
    c = (max(v(:,3))-min(v(:,3)))/2;
end


nv = length(v);
nf = length(f);


%% Disk conformal map for the input surface

% Disk conformal map (using Choi and Lui, J. Sci. Comput. 2015)
disk = disk_conformal_map(v,f);

% (optional) double check the orientation of the disk conformal map 
e1 = [disk(f(:,2),:) - disk(f(:,1),:), disk(f(:,1),1).*0];
e2 = [disk(f(:,3),:) - disk(f(:,1),:), disk(f(:,1),1).*0];
cross12 = cross(e1,e2);
if sum(cross12(:,3) > 0) < sum(cross12(:,3) < 0)
    disp('Corrected orientation.');
    % wrong orientation in flattening
    disk = [disk(:,2), disk(:,1)];
end


%% Search for an optimal Mobius transformation (extending Choi et al., SIAM J. Imaging Sci. 2020)

% Compute the area with normalization
area_v = face_area(f,v); area_v = area_v/sum(area_v);

z = complex(disk(:,1),disk(:,2));

% Fixing the rotation arbitrarily
[~,id_right] = max(v(:,1));
z = z*exp(-1i*angle(z(id_right)));

% Function for calculating the area after the Mobius transformation 
area_map = @(x) face_area(f,stereographic_ellipsoid(...
    x(3)*exp(1i*x(4))*(z-x(1)*exp(1i*x(2)))./(1-conj(x(1)*exp(1i*x(2)))*z),a,b,c))/...
    sum(face_area(f,stereographic_ellipsoid(...
    x(3)*exp(1i*x(4))*(z-x(1)*exp(1i*x(2)))./(1-conj(x(1)*exp(1i*x(2)))*z),a,b,c)));
 

% objective function: 
d_area = @(r) finitemean(abs(log(area_map(r)./area_v)).^2);

% Optimization setup
x0 = [0,0,1,0]; % initial guess, try something diferent if the result is not good
lb = [0,-pi,1e-3,-pi]; % lower bound for the parameters
ub = [1,pi,1e3,pi]; % upper bound for the parameters
options = optimoptions('fmincon','Display','Iter');

% Optimization
x = fmincon(d_area,x0,[],[],[],[],lb,ub,[],options);

% obtain the conformal parameterization with area distortion corrected
fz = x(3)*exp(1i*x(4))*(z-x(1)*exp(1i*x(2)))./(1-conj(x(1)*exp(1i*x(2)))*z);
disk_mobius = [real(fz),imag(fz)];

E_mobius = stereographic_ellipsoid(fz,a,b,c);


%% Construct conformal map from plane to ellispoid using QC composition

% compute Beltrami coefficient
mu = beltrami_coefficient(disk_mobius,f,E_mobius);

% Construct a QC map with the same BC
bd = meshboundaries(f);
disk_map = linear_beltrami_solver(disk_mobius,f,mu,...
    bd{1},disk_mobius(bd{1},:));

% use disk_map to E_mobius to do conformal interpolation
Fx = scatteredInterpolant(disk_map(:,1),disk_map(:,2),E_mobius(:,1));
Fy = scatteredInterpolant(disk_map(:,1),disk_map(:,2),E_mobius(:,2));
Fz = scatteredInterpolant(disk_map(:,1),disk_map(:,2),E_mobius(:,3));


%% Obtain the inverse mapping result
map = zeros(nv,3);
map(:,1) = Fx(disk_mobius(:,1),disk_mobius(:,2));
map(:,2) = Fy(disk_mobius(:,1),disk_mobius(:,2));
map(:,3) = Fz(disk_mobius(:,1),disk_mobius(:,2));

% Finally rotate it to correct the orientation
map = [-map(:,1),-map(:,2),-map(:,3)];


end


function P = stereographic_ellipsoid(p,a,b,c)
% https://mathworld.wolfram.com/Ellipsoid.html
    if size(p, 2) == 1
      p = [real(p), imag(p)];
    end
    u = p(:,1);
    v = p(:,2);
    if size(p,2) < 3
      z = 1 + u.^2 + v.^2;
      P = [2*a*u./z, 2*b*v./z, c*(-1+u.^2+v.^2)./z];
      
      P(isnan(z)|(~isfinite(z)),1) = 0;
      P(isnan(z)|(~isfinite(z)),2) = 0;
      P(isnan(z)|(~isfinite(z)),3) = c;
        
    else
      z = p(:,3);
      P = [(u/a)./(1-z/c), (v/b)./(1-z/c)];
      P(isnan(P)) = Inf;
    end
end


function m = finitemean(A)
    % for avoiding the Inf values caused by division by a very small area
    m = mean(A(isfinite(A)));
end

