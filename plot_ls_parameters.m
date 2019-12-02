function plot_ls_parameters(ls_name, variable_name, save_plot)
% PLOT_LS_PARAMETERS plots landing site parameter map
%
% PLOT_LS_PARAMETERS plots parameters for all fitted parameters and saves
% the outputs
%
% PLOT_LS_PARAMETERS(ls_name) plots parameters for ls_name if ls_name is
% the path to a data file, or plots A0 and H for ls_name is ls_name is the
% name of a simulation environment
%
% PLOT_LS_PARAMETERS(ls_name, variable_name) defines variable to plot
% parameter map of
%
% PLOT_LS_PARAMETERS(ls_name, variable_name, save_plot) controls if plot is
% saved


width = 4;
height = 3;

if nargin < 1 || numel(ls_name) == 0
    
    source_folder = create_static_path('outputs/ls_parameters');
    fprintf('Processing all files in %s\n', source_folder)
    listing = dir(source_folder);
    
    close
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 width height];
    fig.Units = 'inches';
    fig.Position = [fig.Position(1) fig.Position(2) 0 0] + fig.PaperPosition;
    drawnow
    
    for idx = 1:numel(listing)
        if numel(listing(idx).name) > 3 && strcmpi(listing(idx).name(end-3:end), '.mat')
            file_path = sprintf('%s/%s', listing(idx).folder, listing(idx).name);
            fprintf('%s ', listing(idx).name);
            plot_ls_parameters(file_path, [], true);
            fprintf('DONE\n')
        end
    end
    fprintf('ALL DONE\n')
    return
    
end
if nargin < 3 || numel(save_plot) == 0
    save_plot = false;
end

if contains(ls_name, '/') || contains(ls_name, '\')
    load(ls_name, 'data');
    variable_name = data.variable;
    ls_name = strrep(ls_name, '\', '/');
    ls_name = strsplit(ls_name, '/');
    ls_name = ls_name{end};
    ls_name = strrep(ls_name, '_parameter_fit.mat', '');
    ls_name = strrep(ls_name, sprintf('_%s', variable_name), '');
else
    if nargin < 2 || numel(variable_name) == 0
        plot_ls_parameters(ls_name, 'A0', save_plot);
        plot_ls_parameters(ls_name, 'H', save_plot);
        return
    end
    data_path = create_static_path('outputs/ls_parameters');
    data_path = sprintf('%s/%s_%s_parameter_fit.mat', data_path, ls_name, variable_name);
    load(data_path, 'data');
end

clf
plot_h = plot_3d_surface(data.ls_data.ew_matrix/1e3, data.ls_data.ns_matrix/1e3, data.ls_data.elevation_matrix/1e3, data.value_matrix);
plot_h.LineWidth = 0.1;
plot_h.EdgeAlpha = 0.33;
shading faceted

colorbar('off')
c = colorbar('northoutside');
colormap viridis

if strcmp(variable_name, 'H')
    unit_str = ' (m)';
else
    unit_str = '';
end

c.Label.String = sprintf('%s %s fit%s', strrep(ls_name, '_', '\_'), variable_name, unit_str);
c.Position = [0.05, 0.84, 0.9, 0.03];

ax = gca;
ax_position = ax.Position;
ax_position(2) = 0.04;
ax.Position = ax_position;

xh = xlabel('East \rightarrow');
yh = ylabel('\leftarrow North');
zlabel('')
set(xh, 'Rotation', 35);
set(yh, 'Rotation', -35);
ztickformat('%gkm')
xtickformat('%gkm')
ytickformat('%gkm')
zticks([])
view([-45,45])
box on

drawnow

if save_plot
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 width height];
    fig.Units = 'inches';
    fig.Position = [fig.Position(1) fig.Position(2) 0 0] + fig.PaperPosition;
    image_name = create_static_path(sprintf('plots/ls_parameters/%s/%s_%s_fit.png', variable_name, ls_name, variable_name));
    print(image_name, '-dpng', '-r300')
    
    switch variable_name
        case 'H'
            caxis([0,1])
            log_colormap(viridis, 0.5);
        case 'A0'
            caxis([0,0.25])
    end
    drawnow
    image_name = create_static_path(sprintf('plots/ls_parameters/%s/constrained/%s_%s_fit_constrained.png', variable_name, ls_name, variable_name));
    print(image_name, '-dpng', '-r300')
end
end