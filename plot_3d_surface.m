function varargout = plot_3d_surface(x_coords, y_coords, elevation_matrix, color_matrix, title_str, cbar_str, x_label_str, y_label_str, z_label_str, cmap_name)
% PLOT_3D_SURFACE plots 3D surfaces easily - good for representing crater
% surface temperatures etc.
%
% PLOT_3D_SURFACE(x_coords, y_coords, elevation_matrix) plots 3d surface
% with given coordinates and colors defined by elevation
%
% PLOT_3D_SURFACE(x_coords, y_coords, elevation_matrix, color_matrix)
% defines colors of points
%
% PLOT_3D_SURFACE(x_coords, y_coords, elevation_matrix, color_matrix,
% title_str) defines title of plot
%
% PLOT_3D_SURFACE(x_coords, y_coords, elevation_matrix, color_matrix,
% title_str, cbar_str) defines colorbar title string
%
% PLOT_3D_SURFACE(..., x_label_str, y_label_str, z_label_str) defines x, y,
% z axis labels
%
% PLOT_3D_SURFACE(..., x_label_str, y_label_str, z_label_str, cmap_name)
% defines colormap to use
%
% plot_h = PLOT_3D_SURFACE(...) returns handle of plotted surface

if nargin < 5
    title_str = '';
end
if nargin < 6
    cbar_str = '';
end
if nargin < 7
    x_label_str = 'E-W distance (m)';
end
if nargin < 8
    y_label_str = 'N-S distance (m)';
end
if nargin < 9
    z_label_str = 'Elevation (m)';
end
if nargin < 10 || isempty(cmap_name)
    cmap_name = 'viridis';
end
if nargin < 4 || isempty(color_matrix)
    color_matrix = elevation_matrix;
    if nargin < 6
        cbar_str = 'Height';
    end
end

if ~ischar(title_str)
    title_str{2} = strcat('{\rm\color{gray}\fontsize{8}', title_str{2}, '}');
end


clf
colormap(cmap_name);

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 12 6];
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
margin = 0.05;
left = outerpos(1) + ti(1) + margin;
bottom = outerpos(2) + ti(2) + margin;
ax_width = outerpos(3) - ti(1) - ti(3) - 2*margin;
ax_height = outerpos(4) - ti(2) - ti(4) - 2*margin;

if ndims(elevation_matrix) == 3
    num_layers = size(elevation_matrix, 3);
    if ndims(color_matrix) ~= 3
        color_matrix = repmat(color_matrix, 1, 1, num_layers);
    end
    hold on
    for layer_idx = 1:num_layers
        plot_h{layer_idx} = surf(x_coords, y_coords, elevation_matrix(:,:,layer_idx), color_matrix(:,:,layer_idx), 'FaceColor', 'interp');
    end
else
    plot_h = surf(x_coords, y_coords, elevation_matrix, color_matrix, 'FaceColor', 'interp');
end
plot_h.LineWidth = 0.1;
plot_h.EdgeAlpha = 0.33;
axis equal
c = colorbar;
c.Label.String = cbar_str;
xlabel(x_label_str)
ylabel(y_label_str)
zlabel(z_label_str)
title(title_str)

drawnow

if nargout > 0
    varargout{1} = plot_h;
end

end