%*************************************************************************%
% All rights reserved (C) to the authors: Mahmoud SHAQFA, Gary CHOI       %
%                                                                         %
%                                                                         %
% M. Shaqfa Contact:                                                      %
% Department of mechanical engineering,                                   %
% Massachusetts Institute of Technology (MIT)                             %
% Cambridge, MA, USA                                                      %
%               Email: mshaqfa@mit.edu                                    %
% G. Choi Contact:                                                        %
% Department of Mathematics, The Chinese University of Hong Kong,         %
% Hong Kong                                                               %
%               Email: ptchoi@cuhk.edu.hk                                 %
%                                                                         %
%*************************************************************************%
% This code includes implementations for:                                 %
%				- Hemispheroidal harmonics analysis (HSOH)                %
%				- Hemispheroidal parameterization                         %
%                 -- Area-preserving mapping                              %
%                 -- Conformal mapping                                    %
%                 -- Tutte mapping                                        %
%                 -- Balanced mapping                                     %
% This code is part of the paper: "Hemispheroidal parameterization and    %
% harmonic decomposition of simply connected open surfaces"               %
%                                                                         %
%*************************************************************************%
% This library is free software; you can redistribute it and/or modify	  %
% it under the terms of the GNU Lesser General Public License as published%
% by the Free Software Foundation; either version 2.1 of the License, or  %
% (at your option) any later version.                                     %
%                                                                         %
% This library is distributed in the hope that it will be useful,         %
% but WITHOUT ANY WARRANTY; without even the implied warranty of          %
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                    %
% See the GNU Lesser General Public License for more details.        	  %
% You should have received a copy of the GNU Lesser General Public License%
% along with this library; if not, write to the Free Software Foundation, %
% Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA       %
%*************************************************************************%
% Authors of this file: Gary Choi @ 2024


function [bds, outer] = meshboundaries(f)

%  MESHBOUNDARIES Mesh Boundaries
%     [bds, n] = MESHBOUNDARIES(f) extracts boundaries of mesh into cells

nv = max(f(:));
v = zeros(nv, 2);
tr = triangulation(f, v);
fe = tr.freeBoundary;
if isempty(fe)
  bds = {};
  outer = [];
  return;
end
[~, order] = sort(fe(:, 1));
fe = fe(order, :); % reordered boundary edges
ex = zeros(nv, 1); 
vs = false(nv, 1); 
nfe = size(fe, 1);
for i = 1 : nfe
  ex(fe(i, 1)) = fe(i, 2); % table of second edge vertex
end
n = 0;
for i = unique(fe(:))' % for each boundary vertex
  if ~vs(i) 
    n = n + 1; % record boundary components
    [bds{n}, vs] = dfs(ex, vs, i);
  end
end
%% 
m = zeros(n, 1);
for i = 1 : n
  m(i) = size(bds{i}, 1);
end
[~, id] = sort(m, 'descend');
bds = bds(id);
outer = bds{1};

function [bd, vs] = dfs(ex, vs, i)
vs(i) = 1; % marked
bd = i;
if ~vs(ex(i)) % if second vertex hasn't been marked
  [bd, vs] = dfs(ex, vs, ex(i)); % mark second vertex
  bd = [i; bd]; % increase bd
end
