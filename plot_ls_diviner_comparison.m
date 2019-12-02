function plot_h = plot_ls_diviner_comparison(ls_name)
% PLOT_LS_DIVINER_COMPARISON plots comparison summary of landing suite
% Diviner and model temperatures
%
% plot_h = PLOT_LS_DIVINER_COMPARISON(ls_name) plots comparison for given
% environment name or path to model output data file and returns handle for
% plot

global ew_matrix ns_matrix elevation_matrix diviner_data data

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
    load(create_static_path(sprintf('outputs/ls_temperatures/%s_temperatures.mat', ls_name)), 'data');
end

deal_with_jd_overlap = true;

z_arr = data.metadata.z_arr;
lat_arr = data.metadata.crater_data.lat_arr;
long_arr = data.metadata.crater_data.long_arr;
elevation_matrix = data.metadata.crater_data.elevation_matrix;
elevation_matrix = (elevation_matrix-min(elevation_matrix(:)))/1e3;
ew_matrix = data.metadata.crater_data.ew_matrix/1e3;
ns_matrix = data.metadata.crater_data.ns_matrix/1e3;
model_jd_arr = juliandate(data.time_series.dtm_arr);

diviner_ls_name = data.metadata.crater_data.crater_name;
if contains(ls_name, 'compressed')
    diviner_ls_name = sprintf('%s_compressed', diviner_ls_name);
end
diviner_data = load(create_static_path(sprintf('crater_environments/diviner_temperatures/%s_diviner_temperatures.mat', diviner_ls_name)));
diviner_data = diviner_data.output;

avg_error_matrix = NaN(size(elevation_matrix));
abs_error_matrix = NaN(size(elevation_matrix));

for lat_idx = 1:numel(lat_arr)
    lat = lat_arr(lat_idx);
    for long_idx = 1:numel(long_arr)
        long = long_arr(long_idx);
        diviner_bool_arr = diviner_data.lat_arr == lat & diviner_data.long_arr == long & diviner_data.jd_arr >= min(model_jd_arr) & diviner_data.jd_arr <= max(model_jd_arr);
        diviner_T_arr = diviner_data.T_arr(diviner_bool_arr);
        diviner_jd_arr = diviner_data.jd_arr(diviner_bool_arr);
        avg_err_value = 0;
        abs_err_value = 0;
        T_model_arr = squeeze(data.time_series.T_time_series_4dmat(lat_idx, long_idx, 1, :));
        if deal_with_jd_overlap
            if numel(diviner_T_arr) > 0
                overlap_start_T_idx = 1;
                while true
                    jd_value = diviner_jd_arr(overlap_start_T_idx);
                    T_model = interp1(model_jd_arr, T_model_arr, jd_value);
                    overlap_end_T_idx = numel(diviner_T_arr);
                    for T_idx = overlap_start_T_idx:numel(diviner_T_arr)
                        if diviner_jd_arr(T_idx) > jd_value
                            overlap_end_T_idx = T_idx - 1;
                            break
                        end
                    end
                    err_value_arr = T_model - diviner_T_arr(overlap_start_T_idx:overlap_end_T_idx);
                    if max(err_value_arr) > 0 && min(err_value_arr) < 0
                        err_value = 0; % T_model in range of T_diviner
                    else
                        if max(err_value_arr) < 0
                            err_value = max(err_value_arr); % choose err_value closest to 0
                        else
                            err_value = min(err_value_arr);
                        end
                    end
                    avg_err_value = avg_err_value + err_value;
                    abs_err_value = abs_err_value + abs(err_value);
                    
                    overlap_start_T_idx = overlap_end_T_idx + 1;
                    if overlap_start_T_idx > numel(diviner_T_arr)
                        break
                    end
                end
            end
        else
            for T_idx = 1:numel(diviner_T_arr)
                T_model = interp1(model_jd_arr, T_model_arr, diviner_jd_arr(T_idx));
                avg_err_value = avg_err_value + (T_model - diviner_T_arr(T_idx));
                abs_err_value = abs_err_value + abs(T_model - diviner_T_arr(T_idx));
            end
        end
        if deal_with_jd_overlap
            avg_error_matrix(lat_idx, long_idx) = (avg_err_value/numel(unique(diviner_jd_arr)));
            abs_error_matrix(lat_idx, long_idx) = (abs_err_value/numel(unique(diviner_jd_arr)));
        else
            avg_error_matrix(lat_idx, long_idx) = (avg_err_value/numel(diviner_T_arr));
            abs_error_matrix(lat_idx, long_idx) = (abs_err_value/numel(diviner_T_arr));
        end
    end
    
    %     ax_avg = subplot(2,2,1);
    %     plot_map_subroutine(avg_error_matrix);
    %     step_colormap(spectral, 1, 0, true);
    %
    %     ax_abs = subplot(2,2,2);
    %     plot_map_subroutine(abs_error_matrix);
    %     colormap(ax_abs, viridis)
    %     caxis(ax_abs, [0, inf])
    %     drawnow
end


clf

buffer_fraction = 0.25;
lat_idx_buffer = round(numel(lat_arr)*buffer_fraction);
long_idx_buffer = round(numel(long_arr)*buffer_fraction);

error_matrix = abs_error_matrix;

ax_bad = subplot(2,2,3);
max_value = max(max(abs(error_matrix(lat_idx_buffer:end-lat_idx_buffer, long_idx_buffer:end-long_idx_buffer))));
[bad_lat_idx,bad_long_idx] = find(abs(error_matrix) == max_value);
bad_str = sprintf('MAX ERROR | abs. = %0.2fK | avg. = %0.2fK | (lat, long) = (%0.2f, %0.2f)', abs_error_matrix(bad_lat_idx, bad_long_idx), avg_error_matrix(bad_lat_idx, bad_long_idx), lat_arr(bad_lat_idx), long_arr(bad_long_idx));
bad_color = [0.8 0 0];
plot_time_series_subroutine(bad_lat_idx, bad_long_idx, bad_str, bad_color)
x_bad = ew_matrix(bad_lat_idx, bad_long_idx);
y_bad = ns_matrix(bad_lat_idx, bad_long_idx);
z_bad = elevation_matrix(bad_lat_idx,bad_long_idx);

ax_good = subplot(2,2,4);
% max_value = min(min(abs(error_matrix(lat_idx_buffer:end-lat_idx_buffer, long_idx_buffer:end-long_idx_buffer))));
% [good_lat_idx,good_long_idx] = find(abs(error_matrix) == max_value);
good_lat_idx = round(numel(lat_arr)/2 + 0.5);
good_long_idx = round(numel(long_arr)/2 + 0.5);
good_str = sprintf('CENTRE | abs. = %0.2fK | avg. = %0.2fK | (lat, long) = (%0.2f, %0.2f)', abs_error_matrix(good_lat_idx, good_long_idx), avg_error_matrix(good_lat_idx, good_long_idx), lat_arr(good_lat_idx), long_arr(good_long_idx));
good_color = [0 0.6 0];
plot_time_series_subroutine(good_lat_idx, good_long_idx, good_str, good_color)
x_good = ew_matrix(good_lat_idx, good_long_idx);
y_good = ns_matrix(good_lat_idx, good_long_idx);
z_good = elevation_matrix(good_lat_idx,good_long_idx);

ax_abs = subplot(2,2,1);
plot_h = plot_map_subroutine(abs_error_matrix, 'Mean absolute error, T_{Model} - T_{Diviner} (K)');
label_subroutine(true)
colormap(ax_abs, viridis)
caxis(ax_abs, [0, inf])
hold on
scatter3(x_good, y_good, z_good, 100, 'g.')
scatter3(x_bad, y_bad, z_bad, 100, 'r.')


ax_avg = subplot(2,2,2);
plot_map_subroutine(avg_error_matrix, 'Mean error, T_{Model} - T_{Diviner} (K)');
label_subroutine(false)
step_colormap(spectral, 1, 0, true);
hold on
scatter3(x_good, y_good, z_good, 100, good_color, '.')
scatter3(x_bad, y_bad, z_bad, 100, bad_color, '.')


axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');

abs_error_value = mean(mean(abs_error_matrix(lat_idx_buffer:end-lat_idx_buffer, long_idx_buffer:end-long_idx_buffer)));
avg_error_value = mean(mean(avg_error_matrix(lat_idx_buffer:end-lat_idx_buffer, long_idx_buffer:end-long_idx_buffer)));

title_str = sprintf('%s Diviner comparison | abs. = %0.2fK | avg. = %0.2fK', strrep(ls_name, '_', '\_'), abs_error_value, avg_error_value);
text(0.5, 0.98,  title_str, 'FontWeight', 'bold', 'HorizontalAlignment', 'center')



    function plot_h = plot_map_subroutine(plot_matrix, c_label)
        if nargin < 2
            c_label = '';
        end
        plot_h = surf(ew_matrix, ns_matrix, elevation_matrix, plot_matrix);
        plot_h.EdgeAlpha = 0.33;
        plot_h.LineWidth = 0.1;
        axis equal
        box on
        c = colorbar;
        c.Label.String = c_label;
    end

    function plot_time_series_subroutine(lat_idx, long_idx, title_str, color)
        lat = diviner_data.crater_data.lat_arr(lat_idx);
        long = diviner_data.crater_data.long_arr(long_idx);
        
        T_diviner_arr = diviner_data.T_arr(diviner_data.lat_arr == lat & diviner_data.long_arr == long);
        jd_arr = diviner_data.jd_arr(diviner_data.lat_arr == lat & diviner_data.long_arr == long);
        
        T_model_arr = data.time_series.T_time_series_4dmat(...
            data.metadata.crater_data.lat_arr == lat,...
            data.metadata.crater_data.long_arr == long,...
            1, :);
        T_model_arr = squeeze(T_model_arr);
        
        T_model_matrix = data.time_series.T_time_series_4dmat(...
            data.metadata.crater_data.lat_arr == lat,...
            data.metadata.crater_data.long_arr == long,...
            2:end, :);
        T_model_matrix = squeeze(T_model_matrix);
        legend_h_arr = plot(data.time_series.dtm_arr, T_model_matrix, 'Color', 0.85*[1 1 1], 'LineStyle', ':');
        legend_h_arr = legend_h_arr(1);
        hold on
        legend_h_arr(end+1) = plot(data.time_series.dtm_arr, T_model_arr, 'Color', lines(1));
        xlim([data.time_series.dtm_arr(1), data.time_series.dtm_arr(end)])
        legend_h_arr(end+1) = scatter(convert_jd_to_dtm(jd_arr), T_diviner_arr, 'kx');
        % leg = legend(legend_h_arr, {...
        %     'T_{Model} layer     ',...
        %     'T_{Model} surface     ',...
        %     'T_{Diviner}'...
        %     });
        % leg.Orientation = 'horizontal';
        xlabel('Datetime')
        ylabel('Temperature (K)')
        title(title_str, 'Color', color)
        drawnow
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
%         grid off
    end
end