function plot_crater_environments(data_path, variable_name)
% PLOT_CRATER_ENVIRONMENTS plots information about crater environments
%
% PLOT_CRATER_ENVIRONMENTS plots information about all crater environments
% 
% PLOT_CRATER_ENVIRONMENTS(data_path, variable_name) plots information specific crater
% environment for specific variable (e.g. 'Tmax')


if nargin == 0
    source_folder = create_static_path('crater_environments');
    fprintf('Processing all files in %s\n', source_folder)
    listing = dir(source_folder);
    
    width = 4;
    height = 3;
    
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
            
            for variable_name = {'elevation', 'Tmax', 'Tmin'}
                variable_name = variable_name{1};
                plot_crater_environments(file_path, variable_name);
                
                fig = gcf;
                fig.PaperUnits = 'inches';
                fig.PaperPosition = [0 0 width height];
                fig.Units = 'inches';
                fig.Position = [fig.Position(1) fig.Position(2) 0 0] + fig.PaperPosition;
                
                if strcmpi(listing(idx).name(1:2), 'ls')
                    image_name = create_static_path(sprintf('plots/ls_environments/%s/%s', variable_name, strrep(listing(idx).name, '.mat', '.png')));
                else
                    image_name = create_static_path(sprintf('plots/crater_environments/%s/%s', variable_name, strrep(listing(idx).name, '.mat', '.png')));
                end
                print(image_name, '-dpng', '-r300')
            end
            fprintf('DONE\n')
        end
    end
    fprintf('ALL DONE\n')
    return
end

load(data_path, 'data')

crater_name = strsplit(data_path, '/');
crater_name = crater_name{end};
crater_name = strrep(crater_name, '.mat', '');


%% Plotting
clf

ew_matrix = data.ew_matrix/1e3;
ns_matrix = data.ns_matrix/1e3;
elevation_matrix = data.elevation_matrix/1e3;
elevation_matrix = elevation_matrix - min(elevation_matrix(:));

switch variable_name
    case 'elevation'
        value_matrix = elevation_matrix;
        colorbar_string = sprintf('%s elevation (km)', strrep(crater_name, '_', '\_'));
    case 'Tmax'
        if contains(crater_name, 'simulated')
            return
        end
        value_matrix = data.Tmax_matrix;
        colorbar_string = sprintf('%s Diviner Tmax (K)', strrep(crater_name, '_', '\_'));
    case 'Tmin'
        if contains(crater_name, 'simulated')
            return
        end
        value_matrix = data.Tmin_matrix;
        colorbar_string = sprintf('%s Diviner Tmin (K)', strrep(crater_name, '_', '\_'));
end

plot_h = plot_3d_surface(ew_matrix, ns_matrix, elevation_matrix, value_matrix);
plot_h.LineWidth = 0.1;
plot_h.EdgeAlpha = 0.33;

colorbar('off')
c = colorbar('northoutside');
colormap viridis

% hold on
% scatter3(0,0,interp2(elevation_matrix,size(elevation_matrix,2)/2+0.5,size(elevation_matrix,1)/2+0.5), 10, 'r.')

c.Label.String = colorbar_string;
c.Position = [0.05, 0.84, 0.9, 0.03];

ax = gca;
ax_position = ax.Position;
ax_position(2) = 0.04;
ax.Position = ax_position;

xh = xlabel('');
yh = ylabel('\leftarrow North');
zlabel('')
set(xh, 'Rotation', 35);
set(yh, 'Rotation', -35);
ztickformat('%gkm')
xtickformat('%gkm')
ytickformat('%gkm')
xticks(yticks)
yticks(xticks)
zticks([])
view([-45,45])
box on
end