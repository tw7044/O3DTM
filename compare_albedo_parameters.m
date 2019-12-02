% COMPARE_ALBEDO_PARAMETERS plots comparisons of temperatures with
% different albedo parameters

figure(1)
clf
colormap(flipud(viridis))
ax_max = subplot(2,1,1);
plot_compare_model_diviner_temperatures('outputs/1d_simulated_T/simulated_T_seasons_max.mat', false, [0.8500, 0.3250, 0.0980], true);
lim = 40;
y_lim = lim*[-1,1];
text_y = 0.8*lim;
ylim(y_lim)
text(3,text_y, 'A', 'FontWeight', 'bold')
text(8,text_y, 'Maximum temperatures')
ax = gca;
ax.XAxisLocation = 'top';
ax.OuterPosition = [0, 0.5, 1.05, 0.5];
xlabel('Cumulative fraction of lunar surface')
x_tick_labels = {};
for tick = xticks
    x_tick_labels{end+1} = sprintf('%.0f%%', 100*sind(tick));
end
xticklabels(x_tick_labels)
hold on

ax_min = subplot(2,1,2);
plot_compare_model_diviner_temperatures('outputs/1d_simulated_T/simulated_T_seasons_min.mat', false, [0, 0.4470, 0.7410], true);
ylim(y_lim)
text(3,text_y, 'B', 'FontWeight', 'bold')
text(8,text_y, 'Minimum temperatures')
ax = gca;
ax.OuterPosition = [0, 0, 1.05, 0.5];
hold on

figure(2)
clf
ax_ratio = gca;
colormap(flipud(viridis))
[Tmax_matrix, lat_arr, long_arr] = read_diviner_temperatures('max', 4);
Tmin_matrix = read_diviner_temperatures('min', 4);
[long_matrix, lat_matrix] = meshgrid(long_arr, lat_arr);
x_arr = abs(lat_matrix(:));
y_arr = Tmax_matrix(:)./Tmin_matrix(:);
x_bins = min(x_arr):0.25:max(x_arr);
y_bins = (round(min(y_arr))-0.5):0.01:(round(max(y_arr))+0.5);
histogram2(x_arr, y_arr, x_bins, y_bins,'DisplayStyle','tile');
xlabel('Latitude (degrees N/S)')
ylabel('Tmax/Tmin')
title('Ratio of extreme temperatures')
ylim([2,6])

figure(1)
lat_arr = flip(0:2:90);
Tmax_arr = NaN(size(lat_arr));
Tmin_arr = NaN(size(lat_arr));
Tratio_new_arr = NaN(size(lat_arr));
Tratio_old_arr = NaN(size(lat_arr));
hold on
h_max = plot(ax_max, lat_arr, Tmax_arr, 'r', 'LineWidth', 1);
hold on
h_min = plot(ax_min, lat_arr, Tmax_arr, 'r', 'LineWidth', 1);

figure(2)
hold on
h_ratio_new = plot(ax_ratio, lat_arr, Tratio_new_arr, 'r', 'LineWidth', 1);
h_ratio_old = plot(ax_ratio, lat_arr, Tratio_new_arr, 'b', 'LineWidth', 1);
% legend([h_ratio_new, h_ratio_old], {'New parameters', 'Old parameters'})

drawnow

new_parameters = struct;
new_parameters.A0 = 0.08;

old_parameters = struct;
old_parameters.a = 0.06; % coefficinet for calculating albedo
old_parameters.b = 0.25; % coefficient for calculating albedo
old_parameters.A0 = 0.08;

for lat_idx = 1:numel(lat_arr)
    lat = lat_arr(lat_idx);
    
    Tmax_old_arr(lat_idx) = hayne_1d_model(lat, 'seasons', 'max', 0, 0, old_parameters);
    Tmin_old_arr(lat_idx) = hayne_1d_model(lat, 'seasons', 'min', 0, 0, old_parameters);
    
    Tmax_new_arr(lat_idx) = hayne_1d_model(lat, 'seasons', 'max', 0, 0, new_parameters);
    Tmin_new_arr(lat_idx) = hayne_1d_model(lat, 'seasons', 'min', 0, 0, new_parameters);
    
    Tmax_arr(lat_idx) = Tmax_old_arr(lat_idx) - Tmax_new_arr(lat_idx);
    Tmin_arr(lat_idx) = Tmin_old_arr(lat_idx) - Tmin_new_arr(lat_idx);
    Tratio_new_arr(lat_idx) = Tmax_new_arr(lat_idx)/Tmin_new_arr(lat_idx);
    Tratio_old_arr(lat_idx) = Tmax_old_arr(lat_idx)/Tmin_old_arr(lat_idx);
    
    if mod(lat_idx,1) == 0 || lat_idx == numel(lat_arr)
        delete(h_max)
        delete(h_min)
        delete(h_ratio_new)
        delete(h_ratio_old)
        
        h_max = plot(ax_max, lat_arr, Tmax_arr, 'r', 'LineWidth', 1);
        h_min = plot(ax_min, lat_arr, Tmin_arr, 'r', 'LineWidth', 1);
                
        h_ratio_new = plot(ax_ratio, lat_arr, Tratio_new_arr, 'r', 'LineWidth', 1);
        h_ratio_old = plot(ax_ratio, lat_arr, Tratio_old_arr, 'b', 'LineWidth', 1);
        
        legend([h_ratio_new, h_ratio_old], {'New parameters', 'Old parameters'})

        drawnow
    end
end

% figure
% plot(lat_arr, Tmax_arr./Tmin_arr)
% xlabel('Latitude (degrees N/S)')
% ylabel('Tmax change / Tmin change')
% title('Ratio of effect of new albedo parameters on max/min temperatures')
% grid on


