function [map, alpha, beta, gamma, qm_k_i] = opt_map4rec(v,f,map_T,map_C,map_A, max_n, hemispheroid_data)
% A fast method for computing a balanced hemispheroidal map of a simply connected open surface
%
% Input:
% v: nv x 3 vertex coordinates of a simply connected triangle mesh
% f: nf x 3 triangulations of a simply connected triangle mesh
% c: the z-radius of the target hemispheroid (the x,y radii are fixed as 1)
%
% Output:
% map: nv x 3 vertex coordinates of the balanced hemispheroidal map
% alpha,beta,gamma: optimal weighting parameters for Tutte, conformal, and area-preserving
%
% Remark:
% - The input surface should be aligned with the xyz axis beforehand,
%   otherwise the result may be affected
%
% Written by Gary Pui-Tung Choi, 2024

% extract mesh boundary
bd = meshboundaries(f);
if length(bd)~=1
    % if the number of boundaries is not 1, the surface is not simply connected
    error('The input mesh is not a simply connected open surface!');
else
    bdy_v = bd{1};
end

%% compute other mappings and get the Beltrami coefficients
global qm_k_i;% Teh computed coefficients
c = hemispheroid_data.c;

P_map_T = spheroidal_projection(map_T,1,c);
P_map_C = spheroidal_projection(map_C,1,c);
P_map_A = spheroidal_projection(map_A,1,c);
mu_T = zeros(length(f),1); % Tutte
mu_C = beltrami_coefficient(P_map_T,f,P_map_C); % conformal
mu_A = beltrami_coefficient(P_map_T,f,P_map_A); % area-preserving

% do an optimization 
map_w = @(w) run_map(P_map_T,P_map_C,P_map_A,f,bdy_v,mu_T,mu_C,mu_A,w);

% Objective function: optimal reconstruction accuracy
% d = @(w) mean(abs(log(face_area(f,spheroidal_projection(map_w(w),1,c))./face_area(f,v))));
rec_obj = @(w) opt_rec_obj(spheroidal_projection(map_w(w),1,c), v, max_n, hemispheroid_data);

% Optimization setup
x0 = [0.1,0.8,0.1]; % initial guess
lb = [0,0,0]; % lower bound for the parameters
ub = [1,1,1]; % upper bound for the parameters
options = optimoptions('fmincon','Display','iter');

% Optimization (may further supply gradients for better result, not yet implemented)
% x_opt = fmincon(d,x0,[],[],[],[],lb,ub,[],options);
x_opt = fmincon(rec_obj,x0,[],[],[],[],lb,ub,[],options);

alpha = x_opt(1);
beta = x_opt(2);
gamma = x_opt(3);

map_O = map_w(x_opt);
%% Obtain the final hemispheroidal parameterization
map = spheroidal_projection(map_O,1,c);

s = alpha+beta+gamma;
alpha = alpha/s;
beta = beta/s;
gamma = gamma/s;

end

function map_O = run_map(P_map_T,P_map_C,P_map_A,f,bdy_v,mu_T,mu_C,mu_A,w)
    alpha = w(1);
    beta = w(2);
    gamma = w(3);

    % Combine the Beltrami coefficient
    mu_O = alpha*mu_T + beta*mu_C + gamma*mu_A;
    
    % Reconstuct a quasi-conformal map
    if max([alpha,beta,gamma]) == alpha
        map_O = linear_beltrami_solver(P_map_T,f,mu_O,bdy_v,P_map_T(bdy_v,:));
    elseif max([alpha,beta,gamma]) == beta
        map_O = linear_beltrami_solver(P_map_T,f,mu_O,bdy_v,P_map_C(bdy_v,:));
    else 
        map_O = linear_beltrami_solver(P_map_T,f,mu_O,bdy_v,P_map_A(bdy_v,:));
    end
end


function rmse_i = opt_rec_obj(map_i, v, max_n, hemispheroid_data)
    % Objective function to minimize the reconstrcution errors and compute
    % optimal mapping
    % Written by: Mahmoud Shaqfa

    foci = hemispheroid_data.foci;
    zeta = hemispheroid_data.zeta;
    hemispheroid_type = hemispheroid_data.hemispheroid_type;

    % Compute harmonic expansion
    [thetas_i, phis_i] = cart2spheroid(map_i, foci, zeta, hemispheroid_type);
    D_mat_i = hemispheroidal_harmonic_basis(max_n, thetas_i, phis_i, hemispheroid_type, false);
    qm_k_i = D_mat_i\v;
    v_i = real(D_mat_i* qm_k_i);

    % Compute RMSE error
    distances_i = sqrt(sum((v - v_i).^2, 2));
    rmse_i = sqrt(mean(distances_i.^2));
end