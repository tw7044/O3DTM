function plot_h = plot_ls_stability(ls_name, cutoff_mode, cutoff_value)
% PLOT_LS_STABILITY plots stability map for landing site - i.e. if ice will
% be stable below the surface of the regolith
%
% PLOT_LS_STABILITY(ls_name) plots summary of landing site stability for
% given landing site name or data file path
%
% PLOT_LS_STABILITY(ls_name, cutoff_mode) defines which plot variant to use
%
% PLOT_LS_STABILITY(ls_name, cutoff_mode, cutoff_value) defines T value to
% use for cutoff

if nargin < 2 || numel(cutoff_mode) == 0
    cutoff_mode = 'multi';
end
if nargin < 3 || numel(cutoff_value) == 0
    cutoff_value = 112;
end

if contains(ls_name, '/') || contains(ls_name, '\')
    load(ls_name, 'data');
    ls_name = strrep(ls_name, '\', '/');
    ls_name = strsplit(ls_name, '/');
    ls_name = ls_name{end};
    ls_name = strrep(ls_name, 'extremes', '');
    ls_name = strrep(ls_name, '_temperatures', '');
    ls_name = strrep(ls_name, '_.mat', '');
    ls_name = strrep(ls_name, '.mat', '');
else
    load(create_static_path(sprintf('outputs/ls_temperatures_extremes/%s_temperatures_extremes.mat', ls_name)), 'data');
end
T_3dmat = data.extremes.Tmax_3dmat;
z_arr = data.metadata.z_arr;
lat_arr = data.metadata.crater_data.lat_arr;
long_arr = data.metadata.crater_data.long_arr;
elevation_matrix = data.metadata.crater_data.elevation_matrix;
elevation_matrix = (elevation_matrix-min(elevation_matrix(:)))/1e3;
ew_matrix = data.metadata.crater_data.ew_matrix/1e3;
ns_matrix = data.metadata.crater_data.ns_matrix/1e3;

clf

if strcmp(cutoff_mode, 'multi')
    clf
    ax_T = subplot(2,2,1);
    plot_T_min;
    caxis([25, 200])
    step_colormap(spectral, 0.5, cutoff_value);
    label_subroutine(true)
    
    ax_z_top = subplot(2,2,2);
    plot_z(cutoff_value, 'top');
    label_subroutine(false)
    
    ax_z_min = subplot(2,2,3);
    plot_z(cutoff_value, 'min');
    label_subroutine(false)
    
    ax_z_bottom = subplot(2,2,4);
    plot_z(cutoff_value, 'bottom');
    label_subroutine(false)
    
    for ax = [ax_z_top, ax_z_bottom, ax_z_min]
        colormap(ax, flip(viridis))%ramp_colormap(-z_arr, 1))
        caxis(ax, [min(z_arr), max(z_arr)])
        c_handle = ax.Colorbar;
        c_handle.Direction = 'reverse';
    end
%     colormap(ax_z_range, flip(viridis));
%     caxis(ax_z_range, [min(z_arr), max(z_arr)])

    axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    title_str = sprintf('Stability maps for %s with stability temperature = %gK', strrep(ls_name, '_', '\_'), cutoff_value);
    text(0.5, 0.98,  title_str, 'FontWeight', 'bold', 'HorizontalAlignment', 'center')

else
    clf
    colormap viridis
    switch cutoff_mode
        case 'T'
            plot_h = plot_T_min;
            step_colormap(spectral, 0.5, cutoff_value);
            label_subroutine(true)
        case 'z_top'
            plot_h = plot_z(cutoff_value, 'top');
            caxis([min(z_arr), max(z_arr)])
            colormap(flip(viridis))
        case 'z_bottom'
            plot_h = plot_z(cutoff_value, 'bottom');
            caxis([min(z_arr), max(z_arr)])
            colormap(flip(viridis))
        case 'z_range'
            plot_h = plot_z(cutoff_value, 'range');
            caxis([min(z_arr), max(z_arr)])
            colormap(flip(viridis))
        case 'z_min'
            plot_h = plot_z(cutoff_value, 'min');
            caxis([min(z_arr), max(z_arr)])
            colormap(flip(viridis))
    end
    label_subroutine(true);
end


%% SUBROUTINES
    function plot_h = plot_T_min
        plot_matrix = NaN(size(elevation_matrix));
        for lat_idx = 1:numel(lat_arr)
            for long_idx = 1:numel(long_arr)
                plot_matrix(lat_idx, long_idx) = min(T_3dmat(lat_idx, long_idx, :));
            end
        end
        plot_h = plot_subroutine(plot_matrix);
        c = colorbar;
        c.Label.String = sprintf('Minimum sub-surface constant temperature (K)');
    end

    function plot_h = plot_z(cutoff_value, top_bottom_range)
        plot_matrix = NaN(size(elevation_matrix));
        for lat_idx = 1:numel(lat_arr)
            for long_idx = 1:numel(long_arr)
                dz = 0.01; % 1cm
                z_interp_arr = min(z_arr):dz:max(z_arr);
                T_arr = interp1(z_arr, squeeze(T_3dmat(lat_idx, long_idx, :)), z_interp_arr, 'spline');
                z_stable_arr = z_interp_arr(T_arr < cutoff_value);
                if strcmp(top_bottom_range, 'min')
                    [~, value_idx] = min(T_arr);
                    plot_value = z_interp_arr(value_idx);
                elseif numel(z_stable_arr) > 1
                    switch top_bottom_range
                        case 'top'
                            plot_value = min(z_stable_arr);
                        case 'bottom'
                            plot_value = max(z_stable_arr);
                        case 'range'
                            plot_value = max(z_stable_arr) - min(z_stable_arr);
                    end
                else
                    plot_value = NaN;
                end
                plot_matrix(lat_idx, long_idx) = plot_value;
            end
        end
        plot_h = plot_subroutine(plot_matrix);
        c = colorbar;
        switch top_bottom_range
            case 'top'
                c.Label.String = sprintf('Minimum stable depth (m)');
            case 'bottom'
                c.Label.String = sprintf('Maximim stable depth (m)');
            case 'range'
                c.Label.String = sprintf('Total stable depth range (m)');
            case 'min'
                c.Label.String = sprintf('Depth of minimum constant temperature (m)');
        end
    end

    function plot_h = plot_subroutine(plot_matrix)
        plot_h = surf(ew_matrix, ns_matrix, elevation_matrix, plot_matrix, 'FaceColor', 'interp');
        plot_h.LineWidth = 0.1;
        plot_h.EdgeAlpha = 0.33;
        
        axis equal
        shading faceted
    end

    function label_subroutine(axis_labels)
        if axis_labels
%             xlabel('East \rightarrow')
%             ylabel('North \rightarrow')
%             xh = xlabel('East \rightarrow');

            yh = ylabel('\leftarrow North');
%             set(xh, 'Rotation', 35);
            set(yh, 'Rotation', -35);
        end
        zlabel('')
        title('')
        ztickformat('%gkm')
        xtickformat('%gkm')
        ytickformat('%gkm')
        xticks(yticks)
        yticks(xticks)
        zticks([])
        view([-45,45])
%         view([0,90])
        box on
        grid off
    end

end