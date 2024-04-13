function map = hemispheroidal_area_preserving_map(v,f,c)
% Computing hemispheroidal area-preserving map of a simply connected open surface
%
% Input:
% v: nv x 3 vertex coordinates of a simply connected triangle mesh
% f: nf x 3 triangulations of a simply connected triangle mesh
% c: the z-radius of the target hemispheroid (the x,y radii are fixed as 1)
%
% Output:
% map: nv x 3 vertex coordinates of the hemispheroidal area-preserving map
%
% Remark:
% - The input surface should be aligned with the xyz axis beforehand,
%   otherwise the result may be affected
%
% Written by Gary Pui-Tung Choi, 2024

% iteration parameters
epsilon = 0.01;
dt = 0.001;
n_max = 500;

% extract mesh boundary
bd = meshboundaries(f);
if length(bd)~=1
    % if the number of boundaries is not 1, the surface is not simply connected
    error('The input mesh is not a simply connected open surface!');
else
    bdy_v = bd{1};
    bdy_e = [bdy_v, bdy_v([2:end,1])];
end

%% Initial map

% Disk Tutte map with arclength parameterization boundary constraint
bdy_length = sqrt((v(bdy_v,1) - v(bdy_v([2:end,1]),1)).^2 + ...
            (v(bdy_v,2) - v(bdy_v([2:end,1]),2)).^2 + ...
            (v(bdy_v,3) - v(bdy_v([2:end,1]),3)).^2);
partial_edge_sum = zeros(length(bdy_v),1);
for i = 2:length(bdy_length)
    for j = 1:i-1
    partial_edge_sum(i) = partial_edge_sum(i) + bdy_length(j);
    end
end
theta = 2*pi.*partial_edge_sum/sum(bdy_length)+2*pi*0.01;
bdy = exp(theta*1i)';
disk = tutte_map(v,f,bdy_v,bdy); 

% double check the orientation of the disk map
e1 = [disk(f(:,2),:) - disk(f(:,1),:), disk(f(:,1),1).*0];
e2 = [disk(f(:,3),:) - disk(f(:,1),:), disk(f(:,1),1).*0];
cross12 = cross(e1,e2);
if sum(cross12(:,3) > 0) < sum(cross12(:,3) < 0)
    % wrong orientation in flattening
    disk = [disk(:,2), disk(:,1)];
    disp('Corrected orientation.');
end

% Search for an optimal Mobius transformation to further reduce the area distortion 
% (extending the method in Choi et al., SIAM J. Imaging Sci. 2020)

% Compute the face area 
area_v = face_area(f,v); 

z = complex(disk(:,1),disk(:,2));

% Fixing the rotation arbitrarily
[~,id_right] = max(v(:,1));
z = z*exp(-1i*angle(z(id_right)));

% Function for calculating the area after the Mobius transformation 
area_map = @(x) face_area(f,spheroidal_projection(...
    (z-x(1)*exp(1i*x(2)))./(1-conj(x(1)*exp(1i*x(2)))*z),1,c))/...
    sum(face_area(f,spheroidal_projection(...
    (z-x(1)*exp(1i*x(2)))./(1-conj(x(1)*exp(1i*x(2)))*z),1,c)));

% objective function: 
d_area = @(r) finitemean(abs(log(area_map(r)./(area_v/sum(area_v)))).^2);

% Optimization setup
x0 = [0,0]; % initial guess, try something different if the result is not good
lb = [0,-pi]; % lower bound for the parameters
ub = [1,pi]; % upper bound for the parameters
options = optimoptions('fmincon','Display','off');

% Optimization
x = fmincon(d_area,x0,[],[],[],[],lb,ub,[],options);

% obtain the conformal parameterization with area distortion corrected
fz = (z-x(1)*exp(1i*x(2)))./(1-conj(x(1)*exp(1i*x(2)))*z);
map_plane = [real(fz),imag(fz)];

% initial map onto the unit disk
initial_map = map_plane;

%% Hemispheroidal area-preserving map using density diffusion

% define initial density with respect to the hemispheroidal map
map = spheroidal_projection(map_plane,1,c);
rho_f = area_v./face_area(f,map); 
rho_v = f2v_area(map,f)*rho_f; 

step = 1;
rho_v_error = std(rho_v)/mean(rho_v);
disp('Step     std(rho)/mean(rho)');  
disp([num2str(step), '        ',num2str(rho_v_error)]);
    
% perform density equalization via the planar unit disk domain
while rho_v_error > epsilon
    % update rho 
    L = laplace_beltrami(map_plane,f);
    A = lumped_mass_matrix(map_plane,f);
    rho_v_temp = (A+dt*L)\(A*rho_v);
    
    if sum(sum(isnan(rho_v_temp)))~=0
        break;
    end

    % update gradient 
    grad_rho_temp_f = compute_gradient(map_plane,f,rho_v_temp);
    grad_rho_temp_v = f2v_area(map_plane,f)*grad_rho_temp_f;
    
    % projection to ensure that the planar shape remains to be a unit disk
    bdy_normal1 = [-(map_plane(bdy_e(:,1),2) - map_plane(bdy_e(:,2),2)),  (map_plane(bdy_e(:,1),1) - map_plane(bdy_e(:,2),1))];
    bdy_normal1 = [bdy_normal1(:,1) ./ sqrt(bdy_normal1(:,1).^2 + bdy_normal1(:,2).^2) , bdy_normal1(:,2) ./ sqrt(bdy_normal1(:,1).^2 + bdy_normal1(:,2).^2)];
    bdy_normal2 = [-(map_plane(bdy_e([end,1:end-1],1),2) - map_plane(bdy_e([end,1:end-1],2),2)),  (map_plane(bdy_e([end,1:end-1],1),1) - map_plane(bdy_e([end,1:end-1],2),1))];
    bdy_normal2 = [bdy_normal2(:,1) ./ sqrt(bdy_normal2(:,1).^2 + bdy_normal2(:,2).^2) , bdy_normal2(:,2) ./ sqrt(bdy_normal2(:,1).^2 + bdy_normal2(:,2).^2)];
    bdy_normal = (bdy_normal1 + bdy_normal2)/2;

    grad_rho_temp_v(bdy_v,:) = grad_rho_temp_v(bdy_v,:) - ...
        [sum(bdy_normal.*grad_rho_temp_v(bdy_v,:),2).*bdy_normal(:,1), ...
        sum(bdy_normal.*grad_rho_temp_v(bdy_v,:),2).*bdy_normal(:,2)];
    
    % update displacement 
    dmap = -[grad_rho_temp_v(:,1) ./ rho_v_temp, grad_rho_temp_v(:,2) ./ rho_v_temp];
    map_plane = map_plane + dmap*dt;
    map_plane_bdy_norm = sqrt(map_plane(bdy_v,1).^2+map_plane(bdy_v,2).^2);
    map_plane(bdy_v,:) = [map_plane(bdy_v,1)./map_plane_bdy_norm, map_plane(bdy_v,2)./map_plane_bdy_norm];
    
    % ensure bijectivity
    mu = beltrami_coefficient(initial_map,f,map_plane);
    if max(abs(mu))>1
        mu(abs(mu)>1) = mu(abs(mu)>1)./abs(mu(abs(mu)>1))*0.9;
        map_plane = linear_beltrami_solver(initial_map,f,mu,bdy_v,map_plane(bdy_v,:));
    end
    
    step = step + 1;
    rho_v_error_new = std(rho_v_temp)/mean(rho_v_temp);
    disp([num2str(step), '        ',num2str(rho_v_error_new)]);

    % obtain the hemispheroidal map
    map = spheroidal_projection(map_plane,1,c); 

    % re-coupling scheme
    rho_f = area_v./face_area(f,map);
    rho_v = f2v_area(map_plane,f)*rho_f;

    if step > n_max || rho_v_error_new > rho_v_error
        break;
    end
    rho_v_error = rho_v_error_new;
    
end

end


function m = finitemean(A)
    % for avoiding the Inf values caused by division by a very small area
    m = mean(A(isfinite(A)));
end
