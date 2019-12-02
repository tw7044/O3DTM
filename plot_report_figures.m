% PLOT_REPORT_FIGURES plots various figures used in project report
% (uncomment various sections to generate different figures)
 
%% General
clf
colormap(viridis)

%% Density profile
% H_pts = 6;
% H_min = 0.05;
% H_max = 0.30;
% z_min = 0;
% z_max = 0.5;
% z_pts = 100;
% clf
% par = define_parameters;
% rho_d = par.rho_d;
% rho_s = par.rho_s;
% c_arr = viridis(H_pts);
% hold on
% H_arr = linspace(H_min, H_max, H_pts);
% z_arr = linspace(z_min, z_max, z_pts);
% legend_arr = {};
% for H_idx = 1:H_pts
%     H = H_arr(H_idx);
%     rho_arr = [];
%     for z_idx = 1:z_pts
%         z = z_arr(z_idx);
%         rho_arr(z_idx) = rho_d - (rho_d - rho_s)*exp(-z/H);
%     end
%     plot(rho_arr, z_arr*100, 'Color', c_arr(H_idx,:), 'LineWidth', 1)
%     legend_arr{end+1} = sprintf('%.f cm',H*100);
% end
% % c = colorbar;
% % caxis([H_min, H_max])
% % c.Label.String = 'H parameter (m)';
% xlabel('Density (kg m^{-3})')
% ylabel('Depth (cm)')
% set(gca,'Ydir','reverse')
% l = legend(legend_arr, 'Location', 'southwest');
% title(l, 'H parameter')
% xlim([rho_s, rho_d])
% box on
% 
% save_report_figure('density_profile', 3, 4, '-dpng')

%% Temperature variation
% % z_arr = hayne_1d_model(0,0,0,0,0,struct('return_z_grid', true));
% wait_time = 60.5;
% hayne_1d_model(0,wait_time,0,0,0,struct('live_graph', true, 'live_graph_plot_interval', 0));
% ylabel('Layer temperature (K)')
% xlabel('Local time')
% xlim([wait_time-1, wait_time])
% tick_values = linspace(wait_time-1, wait_time, 5);
% xticks(tick_values)
% tick_labels = {};
% for idx = 1:numel(tick_values)
%     tick_value = tick_values(idx);
%     tick_label = sprintf('%02i:00', mod(tick_value+0.5,1)*24);
%     tick_labels{end+1} = tick_label;
% end
% if strcmp(tick_labels{end}, tick_labels{1})
%     tick_labels{end} = '24:00';
% end
% xticklabels(tick_labels);
% box on
% save_report_figure('layer_temperature_variation')

%% Parameter sensitivity
% global t_wait n_pts lat max_range
% lat = 0;
% t_wait = 60;
% n_pts = 100;
% max_range = 0;
% H_arr = linspace(1e-3, 0.30, n_pts);
% emissivity_arr = linspace(0.75, 1, n_pts);
% 
% plot_subroutine('H', H_arr, 'H parameter (cm)')
% save_report_figure('H_sensitivity', 2, 4)
% 
% plot_subroutine('emissivity', emissivity_arr, 'Emissivity')
% save_report_figure('emissivity_sensitivity', 2, 4)
% 
% plot_subroutine('H', H_arr, 'H parameter (cm)')
% save_report_figure('H_sensitivity', 2, 4)

%% Scattering powers
% n_pts = 1e3;
% p_arr = linspace(0,1,5);
% theta_arr = linspace(0, 90, n_pts);
% legend_arr = {};
% colors = viridis(numel(p_arr));
% for p_idx = 1:numel(p_arr)
%     p = p_arr(p_idx);
%     emissivity_arr = (p+1).*cosd(theta_arr).^p.*abs(sind(theta_arr));
%     plot(theta_arr, emissivity_arr, 'Color', colors(p_idx,:), 'LineWidth', 1)
%     legend_arr{end+1} = sprintf('%.2f', p);
%     hold on
% end
% leg = legend(legend_arr, 'Location','best');
% title(leg, 'p')
% xlabel('Angle to surface normal, \theta (degrees)')
% ylabel('Normalised power')
% box on
% xlim([0,90])
% xticks(linspace(0,90,7))
% save_report_figure('scattering_powers')

%% Summary maps
% map_subroutine('outputs/1d_model_maps/H_LOLA-scaled-A0fixed_sloped_60P_0.005Taccuracy_data.mat', 6, 13);
% save_report_figure('H_map', 3.5, 8, '-djpeg')
% 
% map_subroutine('outputs/1d_model_maps/emissivity_LOLA-scaled-A0fixed_sloped_seasons_0.005Taccuracy_data.mat', 0.8, 0.95);
% save_report_figure('emissivity_map', 3.5, 8, '-djpeg')

%% Model - Diviner comparisons
% data_path = 'outputs/1d_simulated_T/simulated_T_60_max.mat';
% plot_compare_model_diviner_temperatures(data_path);
% xlabel('Maximum Diviner temperature (K)')
% ylabel('Maximum 1D model temperature (K)')
% box on
% axis square
% axis_ticks = [100, 200, 300, 400];
% yticks(axis_ticks)
% xticks(axis_ticks)
% annotation_color = 0*ones(1,3);
% annotation('textarrow', [0.3, 0.35], [0.85, 0.7], 'String', 'Shadowing', 'Color', annotation_color)
% annotation('textarrow', [0.55, 0.45], [0.25, 0.4], 'String', 'Scattering', 'Color', annotation_color)
% save_report_figure('1d_diviner_comparison_max', 3.15, 4)
% 
% data_path = 'outputs/1d_simulated_T/simulated_T_seasons_min.mat';
% plot_compare_model_diviner_temperatures(data_path);
% xlabel('Minimum Diviner temperature (K)')
% ylabel('Minimum 1D model temperature (K)')
% box on
% axis square
% % axis_ticks = [100, 200, 300, 400];
% % yticks(axis_ticks)
% % xticks(axis_ticks)
% annotation_color = 0*ones(1,3);
% % annotation('textarrow', [0.3, 0.35], [0.85, 0.7], 'String', 'Shadowing', 'Color', annotation_color)
% % annotation('textarrow', [0.55, 0.45], [0.25, 0.4], 'String', 'Scattering', 'Color', annotation_color)
% save_report_figure('1d_diviner_comparison_min', 3.15, 4)


%% Model - Diviner (lat)
% clf
% colormap(flipud(viridis))
% subplot(2,1,1)
% plot_compare_model_diviner_temperatures('outputs/1d_simulated_T/simulated_T_seasons_max.mat', false, [0.8500, 0.3250, 0.0980], true);
% lim = 40;
% y_lim = lim*[-1,1];
% text_y = 0.8*lim;
% ylim(y_lim)
% text(3,text_y, 'A', 'FontWeight', 'bold')
% text(8,text_y, 'Maximum temperatures')
% ax = gca;
% ax.XAxisLocation = 'top';
% ax.OuterPosition = [0, 0.5, 1.05, 0.5];
% xlabel('Cumulative fraction of lunar surface')
% x_tick_labels = {};
% for tick = xticks
%     x_tick_labels{end+1} = sprintf('%.0f%%', 100*sind(tick));
% end
% xticklabels(x_tick_labels)
% 
% subplot(2,1,2)
% plot_compare_model_diviner_temperatures('outputs/1d_simulated_T/simulated_T_seasons_min.mat', false, [0, 0.4470, 0.7410], true);
% ylim(y_lim)
% text(3,text_y, 'B', 'FontWeight', 'bold')
% text(8,text_y, 'Minimum temperatures')
% ax = gca;
% ax.OuterPosition = [0, 0, 1.05, 0.5];
% save_report_figure('1d_diviner_comparison_lat', 4, 4)
% clf
% axis off
% c = colorbar('north');
% c.Label.String = '';
% c.Position = [0.05, 0.65, 0.9, 0.15];
% c.Ticks = [0, 0.25, 0.5, 0.75 ,1];
% c.TickLabels = {'0', '', 'Number of pixels (normalised values)', '', '1'};
% save_report_figure('1d_diviner_comparison_lat_colorbar', 0.35, 4)

% data = load('outputs/1d_simulated_T/simulated_T_60_max.mat');
% data = data.data;
% Tdiviner_matrix = read_diviner_temperatures(data.T_mode);
% clf
% hold on
% colormap(flipud(viridis))
% marker_size = 1;
% scatter(Tdiviner_matrix(:), data.Tsimulated_matrix(:), marker_size, '.', 'MarkerEdgeAlpha',0.3)
% data = load('outputs/1d_simulated_T/simulated_T_seasons_min.mat');
% data = data.data;
% Tdiviner_matrix = read_diviner_temperatures(data.T_mode);
% scatter(Tdiviner_matrix(:), data.Tsimulated_matrix(:), marker_size, '.', 'MarkerEdgeAlpha',0.3)
% 
% lim = [min(min(xlim), min(ylim)), max(max(xlim), max(ylim))];
% plot(lim, lim, 'k');%'Color', [0.75, 0.75, 0.75])
% xlim(lim)
% ylim(lim)
% box on
% axis square
% save_report_figure('1d_diviner_comparison_both', 3.15, 4)


%% Craters
% simulated_crater_power_comparison_subroutine('simulated_85n', 'max')
% simulated_crater_power_comparison_subroutine('simulated_60n', 'max')
% simulated_crater_power_comparison_subroutine('simulated_30n', 'max')
% simulated_crater_power_comparison_subroutine('simulated_0n', 'max')
% 
% 
% crater_1d_vs_3d_subroutine('outputs/3d/1d_3d_comparison/simulated_30n_10_0.mat', 'crater_1d_3d_comparison');
% crater_1d_vs_3d_subroutine('outputs/3d/1d_3d_comparison/simulated_30n_10_0.74.mat', 'crater_1d_3d_comparison_min');
% c_lim = [-15,0];
% simulated_crater_power_comparison_subroutine('simulated_85n', 'max', c_lim)
% c_lim = caxis;
% crater_p_T_subroutine('simulated_85n', 'max')
% y_lim = ylim;
% save_report_figure('p_T_graph_85n', 2, 4)
% simulated_crater_power_comparison_subroutine('simulated_0n', 'max', c_lim)
% crater_p_T_subroutine('simulated_0n', 'max')
% ylim(y_lim)
% save_report_figure('p_T_graph_0n', 2, 4)
% simulated_crater_power_comparison_subroutine('simulated_30n', 'max', c_lim)
% crater_p_T_subroutine('simulated_30n', 'max')
% ylim(y_lim)
% save_report_figure('p_T_graph_30n', 2, 4)
% simulated_crater_power_comparison_subroutine('simulated_60n', 'max', c_lim)
% crater_p_T_subroutine('simulated_60n', 'max')
% ylim(y_lim)
% save_report_figure('p_T_graph_60n', 2, 4)
% crater_p_T_subroutine('bruce', 'max', 'seasons')




%% Crater 3D
% crater_data = load('outputs/3d/1d_3d_comparison/blagg_match_max.mat');
% crater_data = crater_data.data;
% 
% c_min = min(min(min(crater_data.hayne_1d_3dmat(:,:,1) - crater_data.Tmax_matrix)),...
%     min(min(crater_data.full_3d_3dmat(:,:,1) - crater_data.Tmax_matrix)));
% c_max = max(max(max(crater_data.hayne_1d_3dmat(:,:,1) - crater_data.Tmax_matrix)),...
%     max(max(crater_data.full_3d_3dmat(:,:,1) - crater_data.Tmax_matrix)));
% c_min = -15;
% c_lim = [c_min, c_max];
% 
% local_colors = viridis(256);
% color_extend = round(256*c_max/(-c_min));
% start_color = 1-0.5*(1-local_colors(end, :));
% end_color = [1 1 1];
% start_color = local_colors(end, :);
% for idx = 1:color_extend
%     frac = idx/color_extend;
%     color = start_color - frac*(start_color - end_color);
%     local_colors(end+1, :) = color;
% end
% 
% local_colors = spectral(256);
% cutoff = 128 + round(128*c_max/(-c_min));
% local_colors = local_colors(1:cutoff, :);
% colormap(local_colors)
% 
% plot_matrix = crater_data.hayne_1d_3dmat(:,:,1) - crater_data.Tmax_matrix;
% real_crater_subroutine(crater_data, plot_matrix)
% text(8, 2, '\bfA  \rm1D model error')
% zticks([0,1])
% caxis(c_lim)
% save_report_figure('blagg_error_1d', 2.5)
% 
% plot_matrix = crater_data.full_3d_3dmat(:,:,1) - crater_data.Tmax_matrix;
% real_crater_subroutine(crater_data, plot_matrix)
% text(8, 2, '\bfB  \rm3D model error')
% zticks([0,1])
% caxis(c_lim)
% 
% save_report_figure('blagg_error_3d', 2.5)
% 
% clf
% axis off
% c = colorbar('north');
% caxis(c_lim)
% c.Label.String = 'Model temperature error, T_{Model} - T_{Diviner} (K)';
% c.Position = [0.05, 0.05, 0.9, 0.15];
% save_report_figure('blagg_error_colorbar', 0.6, 4)
% 
% 
% clf
% hold on
% T_1d_matrix = crater_data.hayne_1d_3dmat(:,:,1);
% T_3d_matrix = crater_data.full_3d_3dmat(:,:,1);
% color_arr = flipud(viridis(256));
% diff_arr = linspace(0,max(max(T_3d_matrix - T_1d_matrix)), 256);
% Tdiviner = crater_data.Tmax_matrix;
% for idx = 1:numel(crater_data.Tmax_matrix)
%     Td = Tdiviner(idx);
%     T1d = T_1d_matrix(idx);
%     T3d = T_3d_matrix(idx);
%     color = interp1(diff_arr, color_arr, T3d - T1d);
%     if isnan(color)
%         color = color_arr(1,:);
%     end
%     plot([Td, Td], [T1d, T3d]-Td, 'Color', color)
%     scatter(Td, T3d-Td, 10, color, '^', 'filled')
% end
% text(374.5, 6, '\bfC')
% xlim([0.999*min(Tdiviner(:)), 1.001*max(Tdiviner(:))]);
% plot(xlim, [0,0], 'k--')
% xlabel('Diviner temperature (K)')
% ylabel('T_{Model} - T_{Diviner} (K)')
% box on
% save_report_figure('blagg_error_graph', 2)


% crater_data = load('outputs/3d/1d_3d_comparison/88s_compressed_seasons_max.mat');
% % crater_data = load('../../remote/pc2/outputs/3d/1d_3d_comparison/88s_compressed_seasons_max_CUSTOM_PARAMS.mat');
% crater_data = crater_data.data;
% T_max_matrix = crater_data.full_3d_3dmat(:,:,1);
% crater_data = load('outputs/3d/1d_3d_comparison/88s_compressed_seasons_min.mat');
% % crater_data = load('../../remote/pc2/outputs/3d/1d_3d_comparison/88s_compressed_seasons_min_CUSTOM_PARAMS.mat');
% crater_data = crater_data.data;
% T_min_matrix = crater_data.full_3d_3dmat(:,:,1);
% 
% 
% plot_matrix = T_max_matrix - crater_data.Tmax_matrix;
% c_min = -65; %min(plot_matrix(:));
% c_max = max(plot_matrix(:));
% local_colors = spectral(256);
% cutoff = 128 + round(128*c_max/(-c_min));
% if cutoff <= 256
%     local_colors = local_colors(1:cutoff, :);
% else
%     c_max = -c_min;
% end
% real_crater_subroutine(crater_data, plot_matrix)
% colormap(local_colors)
% caxis([c_min, c_max])
% message = sprintf('\\bfA  \\rmError in model\nmaximim temperature');
% text(9, 1, message)
% zticks([0,1])
% c = colorbar('northoutside');
% c.Label.String = 'Model temperature error, T_{Model} - T_{Diviner} (K)';
% c.Position = [0.05, 0.84, 0.9, 0.03];
% ax = gca;
% ax_position = ax.Position;
% ax_position(2) = 0.04;
% ax.Position = ax_position;
% save_report_figure('88s_error_max', 3.15)
% 
% plot_matrix = T_min_matrix -crater_data.Tmin_matrix;
% c_min = -20; %min(plot_matrix(:));
% c_max = max(plot_matrix(:));
% local_colors = spectral(256);
% cutoff = 128 + round(128*c_max/(-c_min));
% if cutoff <= 256
%     local_colors = local_colors(1:cutoff, :);
% else
%     c_max = -c_min;
% end
% real_crater_subroutine(crater_data, plot_matrix)
% colormap(local_colors)
% caxis([c_min, c_max])
% message = sprintf('\\bfB  \\rmError in model\nminimum temperature');
% text(9, 1, message)
% zticks([0,1])
% c = colorbar('northoutside');
% c.Label.String = 'Model temperature error, T_{Model} - T_{Diviner} (K)';
% c.Position = [0.05, 0.84, 0.9, 0.03];
% ax = gca;
% ax_position = ax.Position;
% ax_position(2) = 0.04;
% ax.Position = ax_position;
% save_report_figure('88s_error_min', 3.15)


% clf
% axis off
% c = colorbar('north');
% caxis(c_lim)
% c.Label.String = 'Model temperature error, T_{Model} - T_{Diviner} (K)';
% c.Position = [0.05, 0.05, 0.9, 0.15];
% save_report_figure('88s_error_colorbar', 0.6, 4)

%% Crater latutitude
% clf
% hold on
% crater_arr = {'simulated_0n', 'simulated_30n', 'simulated_60n', 'simulated_85n'};
% T_mode = 'max';
% lat_idx = 9;
% long_idx = 9;
% legend_arr = {};
% color_arr = flipud(viridis(numel(crater_arr)));
% for crater_idx = 1:numel(crater_arr)
%     crater_name = crater_arr{crater_idx};
%     crater_data = load(sprintf('outputs/3d/simulated_T/%s_%s_seasons_11pts.mat', crater_name, T_mode));
%     crater_data = crater_data.data;
%     x = reshape(crater_data.scattering_power_arr, [], 1);
%     y = reshape(crater_data.Tsimulated_3dmat(:,lat_idx, long_idx), [], 1);
%     T = y(1);
%     y = 100*(y - T)/T;
%     plot(x, y, 'Color', color_arr(crater_idx,:), 'LineWidth', 1)
%     legend_arr{crater_idx} = sprintf('%g%sN', crater_data.center_lat, char(176));
% end
% legend(legend_arr, 'Location', 'southwest');
% ylabel('Relative temperature change')
% ytickformat('%g%%')
% xlabel('Scattering power, p')
% box on
% save_report_figure('p_T_lat_graph_relative')
% 
% clf
% hold on
% for crater_idx = 1:numel(crater_arr)
%     crater_name = crater_arr{crater_idx};
%     crater_data = load(sprintf('outputs/3d/simulated_T/%s_%s_seasons_11pts.mat', crater_name, T_mode));
%     crater_data = crater_data.data;
%     x = reshape(crater_data.scattering_power_arr, [], 1);
%     y = reshape(crater_data.Tsimulated_3dmat(:,lat_idx, long_idx), [], 1);
%     T = y(1);
%     y = (y - T);
%     plot(x, y, 'Color', color_arr(crater_idx,:), 'LineWidth', 1)
%     legend_arr{crater_idx} = sprintf('%g%sN', crater_data.center_lat, char(176));
% end
% legend(legend_arr, 'Location', 'southwest');
% ylabel('Temperature change (K)')
% xlabel('Scattering power, p')
% box on
% save_report_figure('p_T_lat_graph_absolute')



%% Calculate emissivity difference
path = 'outputs/1d_model_maps/emissivity_LOLA-scaled-A0fixed_sloped_seasons_0.005Taccuracy_data.mat';
data = load(path);
data = data.data;
e_matrix = data.value_matrix;
[A0_matrix, lat_arr, long_arr] = read_lola_A0_data;
par = define_parameters;
A0_matrix = scale_data(A0_matrix, par.A0);
A0_cutoff = 0.09;
A0_buffer = 0.01;
e_matrix((abs(lat_arr) > 45) , :) = NaN;

highlands_e_matrix = e_matrix(A0_matrix > A0_cutoff + A0_buffer);
maria_e_matrix = e_matrix(A0_matrix < A0_cutoff - A0_buffer);
highlands_e = nanmean(highlands_e_matrix);
maria_e = nanmean(maria_e_matrix);
difference = maria_e - highlands_e;

maria_std = nanstd(maria_e_matrix);
highlands_std = nanstd(highlands_e_matrix);

fprintf('Highlands:\t%f +/- %f\nMaria:\t\t%f +/- %f\nDifference:\t%f\n', highlands_e, highlands_std, maria_e, maria_std, difference);

plot_matrix = e_matrix;
plot_matrix(A0_matrix > A0_cutoff + A0_buffer) = plot_matrix(A0_matrix >= A0_cutoff + A0_buffer) - 0.5;
plot_matrix(A0_matrix < A0_cutoff - A0_buffer) = plot_matrix(A0_matrix < A0_cutoff - A0_buffer) + 0.5;
plot_map(lat_arr, long_arr, plot_matrix);

%% Subroutines
function crater_p_T_subroutine(crater_name, T_mode, t_wait)
if nargin < 3
    t_wait = 'seasons';
end
crater_data = load(sprintf('outputs/3d/simulated_T/%s_%s_%s_11pts.mat', crater_name, T_mode, t_wait));
crater_data = crater_data.data;
clf
hold on
T_max = max(max(crater_data.Tsimulated_3dmat(1, :, :)));
T_min = min(min(crater_data.Tsimulated_3dmat(1, :, :)));
T_arr = linspace(T_min, T_max, 256);
c_arr = viridis(256);
x = reshape(crater_data.scattering_power_arr, [], 1);
for lat_idx = 1:numel(crater_data.lat_arr)
    for long_idx = 1:numel(crater_data.long_arr)
        y = reshape(crater_data.Tsimulated_3dmat(:,lat_idx, long_idx), [], 1);
        T = y(1);
        y = y - T;
        c_idx = interp1(T_arr, 1:256, T);
        c_idx = round(c_idx);
        c = c_arr(c_idx, :);
        c = [0    0.4470    0.7410 0.25];
%         f = fit(x, y, fittype({'x','x^2', 'x^3'}));
%         plot(x, f.a*x+f.b*x.^2+f.c*x.^3, 'r--')
%         plot(1:3, [f.a, f.b, f.c], 'Color', c)
        plot(x, y, 'Color', c)
    end
end

% c = colorbar;
% caxis([T_min, T_max])
% c.Label.String = 'T_{p=0} (K)';
xticks(0:0.25:1)
ylabel('T_{p} - T_{p=0} (K)')
xlabel('Scattering power, p')
box on
ylim([min(min(crater_data.Tsimulated_3dmat(end,:, :) - crater_data.Tsimulated_3dmat(1,:, :))), max(max(crater_data.Tsimulated_3dmat(end,:, :) - crater_data.Tsimulated_3dmat(1,:, :)))])
end

function crater_1d_vs_3d_subroutine(data_path, output_path)
crater_data = load(data_path);
crater_data = crater_data.data;
plot_matrix = crater_data.full_3d_3dmat - crater_data.hayne_1d_3dmat;
plot_matrix = plot_matrix(:,:,1);
crater_subroutine(crater_data, plot_matrix, 'Temperature difference, T_{3D} - T_{1D} (K)')
save_report_figure(strcat(output_path, '_map'), 3.3)

clf
T_3d_matrix = crater_data.full_3d_3dmat(:,:,1);
T_1d_matrix = crater_data.hayne_1d_3dmat(:,:,1);
% color_matrix = crater_data.ns_matrix;
T_diff_matrix = T_3d_matrix - T_1d_matrix;
scatter(T_1d_matrix(:), T_diff_matrix(:), 100, '.')
ylim([min(T_diff_matrix(:)),max(T_diff_matrix(:))])
xlim([min(T_1d_matrix(:)),max(T_1d_matrix(:))])
box on
xlabel('1D model temperature, T_{1D} (K)')
ylabel('T_{3D} - T_{1D} (K)')
save_report_figure(strcat(output_path, '_graph'), 2)
end

function simulated_crater_power_comparison_subroutine(crater_name, T_mode, c_lim)
crater_data = load(sprintf('outputs/3d/simulated_T/%s_%s_seasons_11pts.mat', crater_name, T_mode));
crater_data = crater_data.data;
idx = find(crater_data.scattering_power_arr == 0.5);
plot_matrix = crater_data.Tsimulated_3dmat(idx,:,:) - crater_data.Tsimulated_3dmat(1,:,:);
crater_subroutine(crater_data, plot_matrix, 'Temperature difference, T_{p=0.5} - T_{p=0} (K)')
if nargin > 2
    caxis(c_lim)
end
save_report_figure(sprintf('crater_%s_%s', crater_name, T_mode), 3.3)
end

function crater_subroutine(crater_data, plot_matrix, colorbar_str)
clf
elevation_matrix = crater_data.elevation_matrix - referenceSphere('moon').Radius;
surf(crater_data.ew_matrix/1e3, crater_data.ns_matrix/1e3, elevation_matrix/1e3, squeeze(plot_matrix))
% zlabel('Elevation (m)')
axis equal
xh = xlabel('\leftarrow East');
yh = ylabel('\leftarrow North');
set(xh, 'Rotation', -35);
set(yh, 'Rotation', 35);
zlim([0,0.5])
zticks([0,0.5])
xticks([-1, 0, 1])
yticks([-1, 0, 1])
ztickformat('%gkm')
xtickformat('%gkm')
ytickformat('%gkm')
shading faceted
view([-135,45])
box on
ax = gca;
ax_position = ax.Position;
ax_position(2) = 0.03;
ax.Position = ax_position;
c = colorbar('northoutside');
c.Label.String = colorbar_str;
c.Position = [0.05, 0.84, 0.9, 0.03];
text(2.25, 0, strcat('\phi =', {' '}, string(crater_data.center_lat), char(176), 'N'))
colormap(viridis)
end

function real_crater_subroutine(crater_data, plot_matrix)
clf
elevation_matrix = crater_data.elevation_matrix - referenceSphere('moon').Radius;
elevation_matrix = elevation_matrix - min(elevation_matrix(:));
surf(crater_data.ew_matrix/1e3, crater_data.ns_matrix/1e3, elevation_matrix/1e3, squeeze(plot_matrix))
% zlabel('Elevation (m)')
axis equal
xh = xlabel('\leftarrow East');
yh = ylabel('\leftarrow North');
set(xh, 'Rotation', -35);
set(yh, 'Rotation', 35);
ztickformat('%gkm')
xtickformat('%gkm')
ytickformat('%gkm')
shading faceted
view([-135,45])
box on
ax = gca;
ax_position = ax.Position;
ax_position(2) = 0.05;
ax_position(4) = 0.9;
ax.Position = ax_position;
end


function map_subroutine(data_path, cmin, cmax)
clf
clim = [cmin,cmax];
data = load(data_path);
data = data.data;
lat_arr = data.lat_arr;
long_arr = data.long_arr;
value_matrix = data.value_matrix;
title_str = data.plot_title;
cbar_str = data.plot_label;
if strcmp(data.variable, 'H')
    value_matrix = value_matrix*100;
    cbar_str = 'H parameter (cm)';
end

label_size = 12;

axes('position', [0.025 0 0.75 1], 'color', 'none', 'box', 'off')
text(-2.5,1,'A','FontSize', label_size)
axesm('mollweid', 'Grid', 'on', 'Frame', 'off');
plot_h = pcolorm(lat_arr, long_arr, value_matrix);
c = colorbar('westoutside');
c.Label.String = cbar_str;
caxis(clim)
box off
axis off
if strcmp(data.variable, 'emissivity')
%     apollo_lat_arr = [0.67408, -3.01239, -3.64530, -8.97301];
%     apollo_long_arr = [23.47297, -23.4215, -17.47136, 15.50019];
%     for mission_idx = 1:numel(apollo_lat_arr)
%         scatterm(apollo_lat_arr(mission_idx), apollo_long_arr(mission_idx), 10, 'ok');
%     end
end


ortho_left = 0.65;
ortho_bottom = 0.035;
ortho_size = 0.45;
ortho_label = 0.8;

axes('position', [ortho_left 1-ortho_bottom-ortho_size ortho_size ortho_size], 'color', 'none', 'box', 'off')
text(ortho_label,ortho_label,'B','FontSize', label_size)
axesm('ortho', 'Grid', 'on', 'Origin', [90,0,0]);
plot_h = pcolorm(lat_arr, long_arr, value_matrix);
caxis(clim)
box off
axis off

axes('position', [ortho_left ortho_bottom ortho_size ortho_size], 'color', 'none', 'box', 'off')
text(ortho_label,-ortho_label,'C','FontSize', label_size)
axesm('ortho', 'Grid', 'on', 'Origin', [-90,0,0]);
plot_h = pcolorm(lat_arr, long_arr, value_matrix);
caxis(clim)
box off
axis off
end



function plot_subroutine(variable_name, variable_arr, x_label)
global t_wait n_pts lat max_range
Tmax_arr = NaN(size(variable_arr));
Tmin_arr = NaN(size(variable_arr));
for idx = 1:n_pts
    Tmax_out = hayne_1d_model(lat, t_wait, 'max', 0, 0, struct(variable_name, variable_arr(idx)));
    Tmax_arr(idx) = Tmax_out(1);
    Tmin_out = hayne_1d_model(lat, t_wait, 'min', 0, 0, struct(variable_name, variable_arr(idx)));
    Tmin_arr(idx) = Tmin_out(1);
end
if strcmp(variable_name, 'H')
    variable_arr = variable_arr*100;
end
clf
hold on
yyaxis left
plot(variable_arr, Tmin_arr, 'LineWidth', 1)
xlabel(x_label)
ylabel('Min. temperature (K)')
if range(ylim) > max_range
    max_range = range(ylim);
end
yyaxis right
plot(variable_arr, Tmax_arr, '--', 'LineWidth', 1)
ylabel('Max. temperature (K)')
if range(ylim) > max_range
    max_range = range(ylim);
end
legend('T_{min}', 'T_{max}')
box on
range_offset = 0.5*[-max_range, max_range];
yyaxis left
ylim(round(mean(ylim)) + range_offset)
yyaxis right
ylim(round(mean(ylim)) + range_offset)
end


