% COMPARE_RESOLUTIONS plots comparison of model values from different
% simulation resolutions

ls_name_1 = 'ls8_16ppd_wide_compressed';
ls_name_2 = 'ls8_128ppd_wide_compressed';
T_mode = 'min';

load(create_static_path(sprintf('outputs/ls_temperatures_extremes/diviner_fit/%s_temperatures_diviner_fit_extremes.mat', ls_name_1)), 'data');
data_1 = data;
load(create_static_path(sprintf('outputs/ls_temperatures_extremes/diviner_fit/%s_temperatures_diviner_fit_extremes.mat', ls_name_2)), 'data');
data_2 = data;

[long_matrix_1, lat_matrix_1] = meshgrid(data_1.metadata.crater_data.long_arr, data_1.metadata.crater_data.lat_arr);
[long_matrix_2, lat_matrix_2] = meshgrid(data_2.metadata.crater_data.long_arr, data_2.metadata.crater_data.lat_arr);

switch T_mode
    case 'max'
        T_matrix_1 = data_1.extremes.Tmax_3dmat(:,:,1);
        T_matrix_2 = data_2.extremes.Tmax_3dmat(:,:,1);
        T_diviner_matrix_1 = data_1.metadata.crater_data.Tmax_matrix;
        T_diviner_matrix_2 = data_2.metadata.crater_data.Tmax_matrix;

    case 'min'
        T_matrix_1 = data_1.extremes.Tmin_3dmat(:,:,1);
        T_matrix_2 = data_2.extremes.Tmin_3dmat(:,:,1);
        T_diviner_matrix_1 = data_1.metadata.crater_data.Tmin_matrix;
        T_diviner_matrix_2 = data_2.metadata.crater_data.Tmin_matrix;
        
    otherwise
        error('Unknown T_mode')
end
T_interp_matrix_1 = interp2(long_matrix_1, lat_matrix_1, T_matrix_1, long_matrix_2, lat_matrix_2);


global ew_matrix ns_matrix elevation_matrix
ew_matrix = data_2.metadata.crater_data.ew_matrix/1e3;
ns_matrix = data_2.metadata.crater_data.ns_matrix/1e3;
elevation_matrix = data_2.metadata.crater_data.elevation_matrix/1e3;
elevation_matrix = elevation_matrix - min(elevation_matrix(:));

figure(1)
clf
ax1 = subplot(2,2,1);
h1 = plot_subroutine(T_interp_matrix_1, 'Interpolated low resolution temperature (K)');

ax2 = subplot(2,2,2);
h2 = plot_subroutine(T_matrix_2, 'High resolution temperature (K)');

ax4 = subplot(2,2,4);
h4 = plot_subroutine(T_interp_matrix_1 - T_matrix_2, 'T_{interpolated low res.} - T_{high res.} (K)');
step_colormap(spectral, 0.5, 0, true, ax4);

ew_matrix = data_1.metadata.crater_data.ew_matrix/1e3;
ns_matrix = data_1.metadata.crater_data.ns_matrix/1e3;
elevation_matrix = data_1.metadata.crater_data.elevation_matrix/1e3;
elevation_matrix = elevation_matrix - min(elevation_matrix(:));

ax3 = subplot(2,2,3);
h3 = plot_subroutine(T_matrix_1, 'Raw low resolution temperature (K)');
xlim(ax3, xlim(ax1))
ylim(ax3, ylim(ax1))
xticks(ax3, xticks(ax1))
yticks(ax3, yticks(ax1))

c_lim = caxis(ax1);
for h = [ax1, ax2, ax3]
    colormap(h, viridis)
    c_lim = [min([c_lim, caxis(h)]), max([c_lim, caxis(h)])];
end
for h = [ax1, ax2, ax3]
    caxis(h, c_lim)
end

axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
title_str = sprintf('Comparison of model resolutions for %s and %s %s temperatures', strrep(ls_name_1, '_', '\_'), strrep(ls_name_2, '_', '\_'), T_mode);
text(0.5, 0.98,  title_str, 'FontWeight', 'bold', 'HorizontalAlignment', 'center')


figure(2)
ew_matrix = data_2.metadata.crater_data.ew_matrix/1e3;
ns_matrix = data_2.metadata.crater_data.ns_matrix/1e3;
elevation_matrix = data_2.metadata.crater_data.elevation_matrix/1e3;
elevation_matrix = elevation_matrix - min(elevation_matrix(:));
clf
ax1 = subplot(2,2,1);
h1 = plot_subroutine(T_diviner_matrix_2, 'T Diviner');

ax2 = subplot(2,2,2);
h2 = plot_subroutine(T_matrix_2 - T_diviner_matrix_2, 'T model error rel. to Diviner');
% step_colormap(spectral, 0.5, 0, true, ax2);
err_value = 50;
divergent_colormap(err_value)



ew_matrix = data_1.metadata.crater_data.ew_matrix/1e3;
ns_matrix = data_1.metadata.crater_data.ns_matrix/1e3;
elevation_matrix = data_1.metadata.crater_data.elevation_matrix/1e3;
elevation_matrix = elevation_matrix - min(elevation_matrix(:));

ax3 = subplot(2,2,3);
h3 = plot_subroutine(T_diviner_matrix_1, 'T Diviner');
xlim(ax3, xlim(ax1))
ylim(ax3, ylim(ax1))
xticks(ax3, xticks(ax1))
yticks(ax3, yticks(ax1))

ax4 = subplot(2,2,4);
h4 = plot_subroutine(T_matrix_1 - T_diviner_matrix_1, 'T model error rel. to Diviner');
% step_colormap(spectral, 0.5, 0, true, ax4);
xlim(ax4, xlim(ax1))
ylim(ax4, ylim(ax1))
xticks(ax4, xticks(ax1))
yticks(ax4, yticks(ax1))

divergent_colormap(err_value)

% c_lim = caxis(ax1);
for h = [ax1, ax3]
    colormap(h, viridis)
%     c_lim = [min([c_lim, caxis(h)]), max([c_lim, caxis(h)])];
end
for h = [ax1, ax3]
    caxis(h, c_lim)
end

axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
title_str = sprintf('DIVINER comparison of model resolutions for %s and %s %s temperatures', strrep(ls_name_1, '_', '\_'), strrep(ls_name_2, '_', '\_'), T_mode);
text(0.5, 0.98,  title_str, 'FontWeight', 'bold', 'HorizontalAlignment', 'center')


drawnow













function h = plot_subroutine(plot_matrix, c_str)
global ew_matrix ns_matrix elevation_matrix
h = surf(ew_matrix, ns_matrix, elevation_matrix, plot_matrix);
h.LineWidth = 0.1;
h.EdgeAlpha = 0.33;
axis equal
% xh = xlabel('East \rightarrow');
yh = ylabel('\leftarrow North');
zlabel('')
title('')
% set(xh, 'Rotation', 35);
set(yh, 'Rotation', -35);
ztickformat('%gkm')
xtickformat('%gkm')
ytickformat('%gkm')
yticks(xticks)
zticks([])
view([-45,45])
box on
c = colorbar;
c.Label.String = c_str;
end