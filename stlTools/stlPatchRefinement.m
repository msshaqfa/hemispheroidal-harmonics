function [vnew, fnew] = stlPatchRefinement(v, f, selected_list)
% Written by Mahmoud Shaqfa, 2021.
% This file to refine a certian selected surface patch
% selected_list: a list of indcies for the selected verts

selected_verts = v(selected_list, :);
selected_faces_list = zeros(1,1);
mid_verts = zeros(1,3);

idx = 1;
for ff = 1:length(f)
    check_mem = ismember(selected_list, f(ff, :));
    if sum(check_mem) == 3
        selected_faces_list(idx) = ff;
        % Calculate mid points per face
        mid_verts(idx, 1:3) = (v(f(ff, 1), :) + v(f(ff, 2), :) + v(f(ff, 3), :))./3;
        idx = idx + 1;
    end
end
mid_verts
[vnew, fnew] = stlAddVerts(v, f, mid_verts);
end