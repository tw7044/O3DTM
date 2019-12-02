function compare_1d_3d_models(crater_name, t_wait, T_mode)
% COMPARE_1D_3D_MODELS saves data file containing lunar crater temperatures
% for 1D and 3D models for a given wait time and temperature mode.
% 
% COMPARE_1D_3D_MODELS(crater_name, t_wait, T_mode) saves data file for
% crater_name with t_wait and T_mode used to define output temperatures,
% e.g. t_wait = 'seasons', T_mode = 'max'

%% Hayne 1d
crater_data =load_crater_environment(crater_name);
crater_data = crater_data.data;
elevation_matrix = crater_data.elevation_matrix;
lat_arr = crater_data.lat_arr;
long_arr = crater_data.long_arr;
ew_matrix = crater_data.ew_matrix;
ns_matrix = crater_data.ns_matrix;
[aspect_matrix, slope_matrix] = convert_height_to_slope(elevation_matrix, lat_arr, long_arr);


hayne_1d_3dmat = NaN([size(elevation_matrix), 21]);
fprintf('Hayne 1d model')
fprintf(repmat('.', 1, numel(lat_arr) - numel('Hayne 1d model') + 4));
fprintf('\n\t\n')
switch T_mode
    case 'max'
        ltim_matrix = crater_data.Tmax_ltim_matrix;
        dtm_matrix = crater_data.Tmax_dtm_matrix;
    case 'min'
        ltim_matrix = crater_data.Tmin_ltim_matrix;
        dtm_matrix = crater_data.Tmin_dtm_matrix;
end
custom_parameters = struct;
if isfield(crater_data, 'A0')
    custom_parameters.A0 = crater_data.A0;
end
for lat_idx =  1:numel(lat_arr)
    lat = lat_arr(lat_idx);
    for long_idx = 1:numel(long_arr)
        aspect = aspect_matrix(lat_idx, long_idx);
        slope = slope_matrix(lat_idx, long_idx);
        switch t_wait
            case 'match'
                t_wait_local = dtm_matrix(lat_idx, long_idx);
                T_mode_local = mod(0.5 + ltim_matrix(lat_idx, long_idx)/24,1);
            case 'match_ltim'
                t_wait_local = 'seasons';
                T_mode_local = T_mode;
            otherwise
                t_wait_local = t_wait;
                T_mode_local = T_mode;
        end
        hayne_1d_3dmat(lat_idx, long_idx, :) = hayne_1d_model(lat, t_wait_local, T_mode_local, aspect, slope, custom_parameters);
    end
    fprintf('\b#\n')
end
%% 3d
full_3d_3dmat = temperature_model_3d(crater_name, T_mode, t_wait, 0, inf);
no_scattering_3d_3dmat = temperature_model_3d(crater_name, T_mode, t_wait, 0, 0);

% elevation_matrix = elevation_matrix - referenceSphere('moon').Radius;
% plot_matrix = full_3d_3dmat - hayne_3dmat;
% plot_matrix = plot_matrix(:,:,1);
% plot_3d_surface(ew_matrix, ns_matrix, elevation_matrix, plot_matrix)
% shading faceted

data = crater_data;
data.hayne_1d_3dmat = hayne_1d_3dmat;
data.full_3d_3dmat = full_3d_3dmat;
data.no_scattering_3d_3dmat = no_scattering_3d_3dmat;
data.t_wait = t_wait;
data.T_mode = T_mode;
data.crater_created = crater_data.created;
data.created = datetime;

data_filename = sprintf('outputs/3d/1d_3d_comparison/%s_%s_%s.mat', crater_name, string(t_wait), string(T_mode));
save(data_filename, 'data')

end