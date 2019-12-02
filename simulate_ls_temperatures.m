function simulate_ls_temperatures(ls_name, parameter)
% SIMULATE_LS_TEMPERATURES simulates landing site temperatures and saves
% results in a data file
%
% SIMULATE_LS_TEMPERATURES(ls_name) specifies landing site to simulate
%
% SIMULATE_LS_TEMPERATURES(ls_name, parameter) specifies parameter to vary

if nargin < 2
    parameter = '';
end
if strcmp(ls_name(1:3), 'ls_')
    for ls_idx = 1:8
        rehash
        simulate_ls_temperatures(sprintf('ls%d%s', ls_idx, ls_name(3:end)), parameter)
    end
    return
end
if numel(parameter) == 0
    data = temperature_model_3d(ls_name);
    target_path = create_static_path(sprintf('outputs/ls_temperatures/%s_temperatures.mat',ls_name));
    save(target_path, 'data')
    clear data
    rehash
    save_ls_diviner_comparison_plot(target_path)
    target_path = remove_time_series_from_file(target_path);
    save_ls_stability_plot(target_path)
    
elseif strcmp(parameter, 'fit') || strcmp(parameter, 'diviner_fit')
    root_filename = create_static_path('outputs/ls_parameters');
    load(sprintf('%s/%s_A0_parameter_fit.mat', root_filename, ls_name), 'data');
    A0_data = data;
    load(sprintf('%s/%s_H_parameter_fit.mat', root_filename, ls_name), 'data');
    H_data = data;
    custom_parameters = struct;
    custom_parameters.A0 = A0_data.value_matrix;
    custom_parameters.H = H_data.value_matrix;
    
    data = temperature_model_3d(ls_name, [], [], [], [], custom_parameters);
    
    data.metadata.fitting_data.A0 = A0_data;
    data.metadata.fitting_data.H = H_data;
    data.metadata.fitting_data.description = 'Modelled temperatures for landing site, using Diviner (max/min) temperatures to constrain model albedo (A0) and H parameter.';
    
    target_path = create_static_path(sprintf('outputs/ls_temperatures/diviner_fit/%s_temperatures_diviner_fit.mat',ls_name));
    save(target_path, 'data')
    clear data
    rehash
    save_ls_diviner_comparison_plot(target_path)
    target_path = remove_time_series_from_file(target_path, 'diviner_fit');
    save_ls_stability_plot(target_path)
else
    switch parameter
        case 'H'
            parameter_value_arr = [0.02, 0.10];  % 0.06 = default parameter
        case 'A0'
            parameter_value_arr = [0.06, 0.12, 0.18];
        case 'Q'
            parameter_value_arr = [0, 0.036]; % 0.018 = default parameter
        case 'old_albedo'
            custom_parameters = struct('a', 0.06, 'b', 0.25);
            data = temperature_model_3d(ls_name, [], [], [], [], custom_parameters);
            rehash
            target_path = create_static_path(sprintf('outputs/ls_temperatures/parameters/%s_temperatures_old_albedo',ls_name));
            save(target_path, 'data')
            clear data
            save_ls_diviner_comparison_plot(target_path)
            target_path = remove_time_series_from_file(target_path, 'parameters');
            save_ls_stability_plot(target_path)
            return
        otherwise
            error('Unexpected parameter value "%s", parameter')
    end
    for parameter_value = parameter_value_arr
        custom_parameters = struct(parameter, parameter_value);
        data = temperature_model_3d(ls_name, [], [], [], [], custom_parameters);
        target_path = create_static_path(sprintf('outputs/ls_temperatures/parameters/%s_temperatures',ls_name));
        target_path = sprintf('%s_%s=%g.mat', target_path, parameter, parameter_value);
        save(target_path, 'data')
        clear data
        rehash
        save_ls_diviner_comparison_plot(target_path)
        target_path = remove_time_series_from_file(target_path, 'parameters');
        save_ls_stability_plot(target_path)
    end
end
end