function varargout = plot_ortho_map(origin_arr, lat_arr, long_arr, value_matrix, title_str, cbar_str, cmap_name)
% PLOT_ORTHO_MAP plots orthographic map of property over lunar surface.
% First argument specifies orientation of orthographic map ('3d' =
% interactive 3d map, -90 = map centered on south pole, [0, 180] = map
% centered on far side etc.). Use second argument to specify path of data
% file (generated by generate_map) to plot map of, or use multiple
% arguments for specific options.
%
% varargout = PLOT_ORTHO_MAP(origin_arr, lat_arr, long_arr, value_matrix,
% title_str, cbar_str, cmap_name)

if nargin < 5
    title_str = '';
end
if nargin < 6
    cbar_str = '';
end
if nargin < 7
    cmap_name = 'viridis';
end
if nargin == 2
    data_path = lat_arr;
    data_from_file = load(data_path);
    data = data_from_file.data;
    lat_arr = data.lat_arr;
    long_arr = data.long_arr;
    value_matrix = data.value_matrix;
    title_str = data.plot_title;
    cbar_str = data.plot_label;
end
if numel(origin_arr) == 0 || strcmp(origin_arr, '3d')
    plot_3d = true;
else
    plot_3d = false;
    if numel(origin_arr) == 1
        origin_arr(2) = 0;
    end
    if numel(origin_arr) == 2
        origin_arr(3) = 0;
    end
    origin_arr = mod(origin_arr + 180, 360) - 180;
end

if ~ischar(title_str)
    title_str{2} = strcat('{\rm\color{gray}\fontsize{8}', title_str{2}, '}');
end
% ensure full wraparound of data
value_matrix(:,end+1) = value_matrix(:,1);
long_arr(end+1) = long_arr(1);

clf
colormap(cmap_name);
if plot_3d
    axesm('globe', 'Grid', 'on', 'ParallelLabel', 'on')
    if sum(sum(isnan(value_matrix))) > 0
        % add white sphere so that NaN values are white rather than pure
        % transparent
        meshm(ones(45,90), [0.25,90,-180], [45,90], -1e-3*ones(45,90), 'FaceColor', 'white');
    end
else
    axesm('ortho', 'Grid', 'on', 'Origin', origin_arr)
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
    ax.Position = [left bottom ax_width ax_height];
end
plot_h = pcolorm(lat_arr, long_arr, value_matrix);
c = contourcbar;
c.Label.String = cbar_str;
title(title_str)

if ~plot_3d
    label_str = sprintf('Rotation: %g%c\nCenter lat: %g%c\nCenter long: %g%c', origin_arr(3), char(176), origin_arr(1), char(176), origin_arr(2), char(176));
    text(-0.98, -0.98, label_str, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom')
    label_r_str = '';
    if origin_arr(1) == 90
        label_r_str = 'North Pole';
    elseif origin_arr(1) == -90
        label_r_str = 'South Pole';
    end
    text(0.98, -0.98, label_r_str, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom')
end
drawnow

if nargout > 0
    [varargout{1:nargout}] = plot_h;
end
end
