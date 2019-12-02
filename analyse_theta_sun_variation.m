% ANALYSE_THETA_SUN_VARIATION plots solar angle values from ray tracing
% code, used for debugging.

load('C:\Users\Oliver\Dropbox\mphys_project\mphys_project_static\crater_environments\ray_tracing\seasons\ls8_16ppd_wide_compressed_ray_tracing.mat'


T_data = load(create_static_path('outputs/ls_temperatures_extremes/ls8_16ppd_wide_compressed_temperatures_extremes.mat'));
T_data = T_data.data;
plot_matrix = T_data.extremes.Tmax_3dmat(:,:,1);


lat_idx = 42;
long_idx = 26;
% plot_matrix(lat_idx, long_idx) = -inf;

theta_arr = 90-rad2deg(squeeze(data.theta_3dmat(lat_idx, long_idx, :)));
fprintf('Fraction of time illuminated for (%i, %i):\t%f\n',lat_idx, long_idx, sum(~isnan(theta_arr))/numel(theta_arr))
fprintf('Diviner extremes:\t %3.1f\t %3.1f\n', T_data.metadata.crater_data.Tmin_matrix(lat_idx,long_idx), T_data.metadata.crater_data.Tmax_matrix(lat_idx,long_idx))
fprintf('Model extremes:  \t %3.1f\t %3.1f\n', T_data.extremes.Tmin_3dmat(lat_idx,long_idx,1), T_data.extremes.Tmax_3dmat(lat_idx,long_idx,1))
fprintf('Model error:     \t%+3.1f\t%+3.1f\n', T_data.extremes.Tmin_3dmat(lat_idx,long_idx,1)-T_data.metadata.crater_data.Tmin_matrix(lat_idx,long_idx), T_data.extremes.Tmax_3dmat(lat_idx,long_idx,1)-T_data.metadata.crater_data.Tmax_matrix(lat_idx,long_idx))

figure(1)
clf
plot_3d_surface(data.ew_matrix, data.ns_matrix, data.elevation_matrix, plot_matrix)
hold on
scatter3(data.ew_matrix(lat_idx,long_idx), data.ns_matrix(lat_idx,long_idx), data.elevation_matrix(lat_idx,long_idx), 100, 'r.')
drawnow


figure(2)
clf
hold on
plot(data.dtm_arr, theta_arr)
no_sun_arr = NaN(size(theta_arr));
no_sun_arr(isnan(theta_arr)) = 0;
plot(data.dtm_arr, no_sun_arr, 'LineWidth',2)
xlabel('DTM')
ylabel('Angle of sun (0 = on *surface* local horizon, 90 = perpendicular to surface)') 
drawnow

