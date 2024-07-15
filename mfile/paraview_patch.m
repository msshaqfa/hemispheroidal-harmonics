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
% Authors of this file: Mahmoud S. Shaqfa @ 2024



function fig = paraview_patch(v, f, maps)
% Plot STL meshes with Paraview maps. Mahmoud Shaqfa (EPFL)

fig = figure; 
% set(fig,'renderer','painters'); % To save vectorized SVG and PDF formats

if nargin < 3
    patch('Faces',f,'Vertices',v,'FaceColor',[42,63,188]./255,'LineWidth',0.5);
else
%     patch('Faces',f,'Vertices',v,'FaceColor','interp','FaceVertexCData',maps,...
%         'EdgeColor','k', 'LineWidth', 0.1, 'LineStyle', '-');
    patch('Faces',f,'Vertices',v,'FaceColor','interp','FaceVertexCData',maps,...
        'EdgeColor','k', 'LineStyle', 'none');
%     
    % Define a color map similar to Paraview's shades
    red_color = zeros(1, 255); green_color = zeros(1, 255); blue_color = zeros(1, 255);
    red_color(1:127) = linspace(42,220,127); red_color(128:255) = linspace(220,174,127+1);
    green_color(1:127) = linspace(63,220,127); green_color(128:255) = linspace(220,0,127+1);
    blue_color(1:127) = linspace(181,220,127); blue_color(128:255) = linspace(220,22,127+1);
    
    ParaviewMap = [red_color', green_color', blue_color']./255;
    
    colormap(ParaviewMap);
%     colormap(winter);
%     colormap default;
    set(gcf,'color','w');
    cb = colorbar;
    set(cb,'position',[.15 .1 .05 .2])
end
axis equal tight off
view(45, 45);
gca.Clipping = 'off';
