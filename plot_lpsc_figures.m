clf
colormap viridis

ls_name = 'ls8_16ppd_wide_compressed';
load(create_static_path(sprintf('outputs/ls_temperatures/diviner_fit/%s_temperatures_diviner_fit.mat', ls_name)), 'data');
diviner_ls_name = data.metadata.crater_data.crater_name;
if contains(ls_name, 'compressed')
    diviner_ls_name = sprintf('%s_compressed', diviner_ls_name);
end
diviner_data = load(create_static_path(sprintf('crater_environments/diviner_temperatures/%s_diviner_temperatures.mat', diviner_ls_name)));
diviner_data = diviner_data.output;
z_arr = data.metadata.z_arr;
lat_arr = data.metadata.crater_data.lat_arr;
long_arr = data.metadata.crater_data.long_arr;
elevation_matrix = data.metadata.crater_data.elevation_matrix;
elevation_matrix = (elevation_matrix-min(elevation_matrix(:)));
ew_matrix = data.metadata.crater_data.ew_matrix;
ns_matrix = data.metadata.crater_data.ns_matrix;
model_jd_arr = juliandate(data.time_series.dtm_arr);

deal_with_jd_overlap = true;

lat_idx = 33;
long_idx = 21;
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


[aspect_matrix, slope_matrix] = convert_height_to_slope(data.metadata.crater_data.elevation_matrix,...
                                                        data.metadata.crater_data.lat_arr,...
                                                        data.metadata.crater_data.long_arr);


custom_params = struct;
d_H = load(create_static_path('outputs/ls_parameters/ls8_16ppd_wide_compressed_H_parameter_fit.mat'));
d_A0 = load(create_static_path('outputs/ls_parameters/ls8_16ppd_wide_compressed_A0_parameter_fit.mat'));
custom_params.H = d_H.data.value_matrix(lat_idx, long_idx);
custom_params.A0 = d_A0.data.value_matrix(lat_idx, long_idx);


data_1d = hayne_1d_model(lat, 'seasons', 'times', aspect_matrix(lat_idx,long_idx), slope_matrix(lat_idx,long_idx), custom_params, [], long);

clf
hold on
line_color = lines;
line_color = line_color(1,:);
plot(data.time_series.dtm_arr, T_model_arr, 'Color', line_color);
plot(data_1d.dtm_arr, data_1d.T_arr, '-');
xlim([data.time_series.dtm_arr(1), data.time_series.dtm_arr(end)])
scatter(convert_jd_to_dtm(jd_arr), T_diviner_arr, 50, 'k.');
box on
xlabel('Date')
ylabel('Surface temperature (K)')
save_paper_figure('lpsc/ls_diviner_comparison_time_series')
% 
% avg_error_matrix = NaN(size(elevation_matrix));
% abs_error_matrix = NaN(size(elevation_matrix));
% for lat_idx = 1:numel(lat_arr)
%     lat = lat_arr(lat_idx);
%     for long_idx = 1:numel(long_arr)
%         long = long_arr(long_idx);
%         diviner_bool_arr = diviner_data.lat_arr == lat & diviner_data.long_arr == long & diviner_data.jd_arr >= min(model_jd_arr) & diviner_data.jd_arr <= max(model_jd_arr);
%         diviner_T_arr = diviner_data.T_arr(diviner_bool_arr);
%         diviner_jd_arr = diviner_data.jd_arr(diviner_bool_arr);
%         avg_err_value = 0;
%         abs_err_value = 0;
%         T_model_arr = squeeze(data.time_series.T_time_series_4dmat(lat_idx, long_idx, 1, :));
%         if deal_with_jd_overlap
%             if numel(diviner_T_arr) > 0
%                 overlap_start_T_idx = 1;
%                 while true
%                     jd_value = diviner_jd_arr(overlap_start_T_idx);
%                     T_model = interp1(model_jd_arr, T_model_arr, jd_value);
%                     overlap_end_T_idx = numel(diviner_T_arr);
%                     for T_idx = overlap_start_T_idx:numel(diviner_T_arr)
%                         if diviner_jd_arr(T_idx) > jd_value
%                             overlap_end_T_idx = T_idx - 1;
%                             break
%                         end
%                     end
%                     err_value_arr = T_model - diviner_T_arr(overlap_start_T_idx:overlap_end_T_idx);
%                     if max(err_value_arr) > 0 && min(err_value_arr) < 0
%                         err_value = 0; % T_model in range of T_diviner
%                     else
%                         if max(err_value_arr) < 0
%                             err_value = max(err_value_arr); % choose err_value closest to 0
%                         else
%                             err_value = min(err_value_arr);
%                         end
%                     end
%                     avg_err_value = avg_err_value + err_value;
%                     abs_err_value = abs_err_value + abs(err_value);
%                     
%                     overlap_start_T_idx = overlap_end_T_idx + 1;
%                     if overlap_start_T_idx > numel(diviner_T_arr)
%                         break
%                     end
%                 end
%             end
%         else
%             for T_idx = 1:numel(diviner_T_arr)
%                 T_model = interp1(model_jd_arr, T_model_arr, diviner_jd_arr(T_idx));
%                 avg_err_value = avg_err_value + (T_model - diviner_T_arr(T_idx));
%                 abs_err_value = abs_err_value + abs(T_model - diviner_T_arr(T_idx));
%             end
%         end
%         if deal_with_jd_overlap
%             avg_error_matrix(lat_idx, long_idx) = (avg_err_value/numel(unique(diviner_jd_arr)));
%             abs_error_matrix(lat_idx, long_idx) = (abs_err_value/numel(unique(diviner_jd_arr)));
%         else
%             avg_error_matrix(lat_idx, long_idx) = (avg_err_value/numel(diviner_T_arr));
%             abs_error_matrix(lat_idx, long_idx) = (abs_err_value/numel(diviner_T_arr));
%         end
%     end
% end
% 
% clf
% plot_subroutine(ew_matrix, ns_matrix, elevation_matrix, abs_error_matrix);
% label_subroutine
% ls_loc_subroutine
% c = colorbar('northoutside');
% c.Label.String = sprintf('Absolute model temperature deviation (K)');
% c.Position = [0.05, 0.84, 0.9, 0.03];
% ax = gca;
% ax_position = ax.Position;
% ax_position(2) = 0.03;
% ax.Position = ax_position;
% save_paper_figure('lpsc/ls_diviner_comparison_abs_error_full_colorbar')
% caxis([0,30])
% c.TickLabels(end) = {'>30'};
% save_paper_figure('lpsc/ls_diviner_comparison_abs_error_tight_colorbar')
% 
% 
% 
% 
% 
% 
% function label_subroutine
% % xh = xlabel('East \rightarrow');
% xticks('auto')
% yticks('auto')
% yh = ylabel('\leftarrow North');
% zlabel('')
% xlabel('')
% title('')
% % set(xh, 'Rotation', 35);
% set(yh, 'Rotation', -35);
% ztickformat('%gkm')
% xtickformat('%gkm')
% ytickformat('%gkm')
% yticks(xticks)
% if numel(zticks) == 0
%     zticks([0 10])
% end
% zticks([min(zticks), max(zticks)])
% view([-45,45])
% box on
% end
% 
% 
% function prepare_paper_size(height, width)
% if nargin < 1 || isempty(height)
%     height = 3;
% end
% if nargin < 2 || isempty(width)
%     width = 4;
% end
% fig = gcf;
% fig.PaperUnits = 'inches';
% fig.PaperPosition = [0 0 width height];
% fig.Units = 'inches';
% fig.Position = [fig.Position(1) fig.Position(2) 0 0] + fig.PaperPosition;
% drawnow
% end
% 
% function parameter_sensitivity_subroutine(variable_name, variable_arr, x_label)
% global t_wait n_pts lat max_range
% Tmax_arr = NaN(size(variable_arr));
% Tmin_arr = NaN(size(variable_arr));
% for idx = 1:round(n_pts/2)
%     fprintf('_')
% end
% fprintf('\n')
% for idx = 1:n_pts
%     if mod(idx, 2) == 0
%         fprintf('\b:')
%     else
%         fprintf('.')
%     end
%     Tmax_out = hayne_1d_model(lat, t_wait, 'max', 0, 0, struct(variable_name, variable_arr(idx)));
%     Tmax_arr(idx) = Tmax_out(1);
%     Tmin_out = hayne_1d_model(lat, t_wait, 'min', 0, 0, struct(variable_name, variable_arr(idx)));
%     Tmin_arr(idx) = Tmin_out(1);
% end
% fprintf('\n')
% if strcmp(variable_name, 'H')
%     variable_arr = variable_arr*100;
% end
% clf
% hold on
% yyaxis left
% plot(variable_arr, Tmin_arr, 'LineWidth', 1)
% xlabel(x_label)
% ylabel('Min. temperature (K)')
% if range(Tmin_arr) > max_range
%     max_range = ceil(range(Tmin_arr)+1);
% end
% yyaxis right
% plot(variable_arr, Tmax_arr, '--', 'LineWidth', 1)
% ylabel('Max. temperature (K)')
% if range(Tmax_arr) > max_range
%     max_range = ceil(range(Tmax_arr)+1);
% end
% legend('T_{min}', 'T_{max}')
% box on
% range_offset = 0.5*[-max_range, max_range];
% yyaxis left
% ylim(round(min(Tmin_arr) + 0.5*range(Tmin_arr)) + range_offset)
% yyaxis right
% ylim(round(min(Tmax_arr) + 0.5*range(Tmax_arr)) + range_offset)
% end
% 
% function real_crater_subroutine(crater_data, plot_matrix)
% clf
% elevation_matrix = crater_data.elevation_matrix - referenceSphere('moon').Radius;
% elevation_matrix = elevation_matrix - min(elevation_matrix(:));
% plot_3d_surface(crater_data.ew_matrix/1e3, crater_data.ns_matrix/1e3, elevation_matrix/1e3, squeeze(plot_matrix))
% colorbar('off')
% grid off
% label_subroutine
% shading faceted
% box on
% ax = gca;
% ax_position = ax.Position;
% ax_position(2) = 0.05;
% ax_position(4) = 0.9;
% ax.Position = ax_position;
% end
% 
% 
% function ls_loc_subroutine(style, ls_name)
% if nargin == 0
%     style = 'r.';
% end
% if nargin < 2
%     ls_name = 'ls8_16ppd_wide_compressed';
% end
% data = load_crater_environment(ls_name);
% ls_lat = -82.7;
% ls_long = 33.50;
% lat_idx = interp1(data.lat_arr, 1:numel(data.lat_arr), ls_lat);
% long_idx = interp1(data.long_arr, 1:numel(data.long_arr), ls_long);
% % fprintf('LAT  =%g\t min=%g\t max=%g\t len=%g\t idx=%g\n', ls_lat, min(lat_arr), max(lat_arr), numel(lat_arr), lat_idx)
% % fprintf('LONG =%g\t min=%g\t max=%g\t len=%g\t idx=%g\n', ls_long, min(long_arr), max(long_arr), numel(long_arr), long_idx)
% 
% elevation_matrix = data.elevation_matrix - min(data.elevation_matrix(:));
% elev = interp2(elevation_matrix, long_idx, lat_idx);
% ew = interp2(data.ew_matrix, long_idx, lat_idx);
% ns = interp2(data.ns_matrix, long_idx, lat_idx);
% 
% elev = elev/1e3;
% ew = ew/1e3;
% ns = ns/1e3;
% hold on
% scatter3(ew, ns, elev, 25, style)
% hold off
% end
% 
% function plot_h = plot_subroutine(ew_matrix, ns_matrix, elevation_matrix, plot_matrix)
% elevation_matrix = elevation_matrix - min(elevation_matrix(:));
% plot_h = surf(ew_matrix/1e3, ns_matrix/1e3, elevation_matrix/1e3, plot_matrix, 'FaceColor', 'interp');
% plot_h.LineWidth = 0.1;
% plot_h.EdgeAlpha = 0.33;
% grid off
% axis equal
% shading faceted
% label_subroutine
% end
% 
