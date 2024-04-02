function A = lumped_mass_matrix(v,f)
% compute the lumped mass matrix for FEM Laplacian
% A(i,j) = \int \phi_i \phi_j 
% =    1-ring area / 6 if i = j
%   (|T_1| + |T_2| /12 if i ~= j

nv = length(v);
f1 = f(:,1); f2 = f(:,2); f3 = f(:,3);

% edge length
l = [sqrt(sum((v(f2,:) - v(f3,:)).^2,2)), sqrt(sum((v(f3,:) - v(f1,:)).^2,2)), sqrt(sum((v(f1,:) - v(f2,:)).^2,2))];
l1 = l(:,1); l2 = l(:,2); l3 = l(:,3);

% Heron's formula
s = (l1 + l2 + l3)*0.5;
area = sqrt( s.*(s-l1).*(s-l2).*(s-l3));
 
% construct matrix
II = [f1; f2; f3];
JJ = [f1; f2; f3];
V = [area; area; area]/3;
A = sparse(II,JJ,V,nv,nv);


end