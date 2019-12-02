function [data_filename, data] = fit_ls_parameters(ls_name, variable_name, fixed_data_paths, break_idx, run_model)
% FIT_LS_PARAMETERS uses max/min diviner temperature data to find
% parameters values for a specific landing site
%
% [data_filename, data] = FIT_LS_PARAMETERS(...) saves parameter values and
% returns filename of saved data file and copy of data file
%
% ... = FIT_LS_PARAMETERS(ls_name) fits parameters for environment with
% name ls_name for A0 then H
%
% ... = FIT_LS_PARAMETERS(ls_name, variable_name) fits parameters for
% environment for parameter variabe_name
% 
% ... = FIT_LS_PARAMETERS(ls_name, variable_name, fixed_data_paths) fits
% parameters using parameters already fitted with data files stored
% fixed_data_paths
% 
% ... = FIT_LS_PARAMETERS(ls_name, variable_name, fixed_data_paths,
% break_idx) fits parameters with number of iterations defined by break_idx
% (default = 5)
% 
% ... = FIT_LS_PARAMETERS(ls_name, variable_name, fixed_data_paths,
% break_idx, run_model) defines if model should be run at end to simulate
% temperatures with final parameter values



if strcmp(ls_name(1:3), 'ls_')
    for ls_idx = 1:8
        rehash
        switch nargin
            case 1
                fit_ls_parameters(sprintf('ls%d%s', ls_idx, ls_name(3:end)));
            case 2
                fit_ls_parameters(sprintf('ls%d%s', ls_idx, ls_name(3:end)), variable_name);
            case 3
                fit_ls_parameters(sprintf('ls%d%s', ls_idx, ls_name(3:end)), variable_name, fixed_data_paths);
            case 4
                fit_ls_parameters(sprintf('ls%d%s', ls_idx, ls_name(3:end)), variable_name, fixed_data_paths, break_idx);
            otherwise
                warning('Unknown inputs')
        end
    end
    return
end
if nargin < 2 || numel(variable_name) == 0
    fprintf('\n\n== Running A0 & H fitting for %s ==\n\n', ls_name)
    fprintf('Fitting for A0\n')
    A0_filename = fit_ls_parameters(ls_name, 'A0', [], [], false);
    fprintf('\nFitting for H\n')
    fit_ls_parameters(ls_name, 'H', A0_filename, [], true);
    return
end
if nargin < 3 || numel(fixed_data_paths) == 0
    fixed_data_paths = {};
end
if ischar(fixed_data_paths)
    fixed_data_paths = {fixed_data_paths};
end
if nargin < 4 || numel(break_idx) == 0
    break_idx = 5;
end
if nargin < 5 || numel(run_model) == 0
    run_model = true;
end

parameters = define_parameters();

t_wait = 'seasons';
debug = false;

switch variable_name
    case 'H'
        value_initial_arr = [0, parameters.H, 0.50,];
        value_limits = [0, inf];
        T_mode = 'min';
    case 'A0'
        value_initial_arr = [0, parameters.A0, 0.5];
        value_limits = [0, 0.999*(1 - parameters.a*2^3 - parameters.b)];
        T_mode = 'max';
    otherwise
        error('Unknown variable chosen')
end

ls_data = load_crater_environment(ls_name);

switch T_mode
    case 'max'
        T_diviner_matrix = ls_data.Tmax_matrix;
    case 'min'
        T_diviner_matrix = ls_data.Tmin_matrix;
end

T_model_3dmat = NaN(size(ls_data.elevation_matrix));
value_3dmat = value_initial_arr(2)*ones(size(ls_data.elevation_matrix));
for loop_idx = 1:break_idx
    fprintf('Model run %i: ', loop_idx)
    if debug && loop_idx == 1
        figure(1)
        plot_3d_surface(ls_data.ew_matrix, ls_data.ns_matrix, ls_data.elevation_matrix, value_3dmat(:,:,end), 'Current values (before run #1)')
        figure(2)
        plot_3d_surface(ls_data.ew_matrix, ls_data.ns_matrix, ls_data.elevation_matrix, NaN(size(ls_data.elevation_matrix)), 'Previous values')
        figure(3)
        plot_3d_surface(ls_data.ew_matrix, ls_data.ns_matrix, ls_data.elevation_matrix, NaN(size(ls_data.elevation_matrix)), 'Value change')
        figure(4)
        plot_3d_surface(ls_data.ew_matrix, ls_data.ns_matrix, ls_data.elevation_matrix, NaN(size(ls_data.elevation_matrix)), 'Current T')
        figure(5)
        plot_3d_surface(ls_data.ew_matrix, ls_data.ns_matrix, ls_data.elevation_matrix, NaN(size(ls_data.elevation_matrix)), 'Current T - Diviner T')
        drawnow
    end
    
    custom_parameters = struct(variable_name, value_3dmat(:,:,loop_idx));
    for parameter_idx = 1:numel(fixed_data_paths)
        fixed_data = load(fixed_data_paths{parameter_idx});
        fixed_data = fixed_data.data;
        custom_parameters.(fixed_data.variable) = fixed_data.value_matrix;
    end
    T_model_output = temperature_model_3d(ls_name, T_mode, t_wait, [], [], custom_parameters);
    T_model_3dmat(:,:,loop_idx) = T_model_output(:,:,1); % only keep surface temperatures
    
    % Update values
    for lat_idx = 1:numel(ls_data.lat_arr)
        for long_idx = 1:numel(ls_data.long_arr)
            T_diviner = T_diviner_matrix(lat_idx, long_idx);
            T_model = T_model_3dmat(lat_idx, long_idx, loop_idx);
            if loop_idx == 1
                if T_diviner == T_model
                    value_3dmat(lat_idx, long_idx, loop_idx+1) = value_initial_arr(2);
                elseif T_diviner > T_model
                    value_3dmat(lat_idx, long_idx, loop_idx+1) = value_initial_arr(1);
                else
                    value_3dmat(lat_idx, long_idx, loop_idx+1) = value_initial_arr(3);
                end
            else
                T_model_interp_arr = squeeze(T_model_3dmat(lat_idx, long_idx, loop_idx-1:loop_idx));
                value_interp_arr = squeeze(value_3dmat(lat_idx, long_idx, loop_idx-1:loop_idx));
                if numel(unique(value_interp_arr)) == 1
                    new_value = value_interp_arr(1); % deal with constant value
                    if debug
                        fprintf('unique value\t%i\t%i\t%g > %g\n', lat_idx, long_idx, value_interp_arr(2), new_value)
                    end
                elseif numel(unique(T_model_interp_arr)) == 1
                    if T_diviner == T_model
                        new_value = value_interp_arr(2);
                    elseif T_diviner > T_model
                        new_value = value_initial_arr(1);
                    else
                        new_value = value_initial_arr(3);
                    end
                    if debug
                        fprintf('unique T\t\t%i\t%i\t%g > %g\n', lat_idx, long_idx, value_interp_arr(2), new_value)
                    end
                elseif sign(T_model_interp_arr(2) - T_model_interp_arr(1))*sign(value_interp_arr(2) - value_interp_arr(1)) > 0
                    % deal with values which do not follow expected trend
                    % (expect value increase to decrease T) due to
                    % environment changing etc.
                    new_value = value_interp_arr(2);
                    if debug
                        fprintf('sign\t\t%i\t%i\t%g > %g\n', lat_idx, long_idx, value_interp_arr(2), new_value)
                        xx = interp1(T_model_interp_arr, value_interp_arr, T_diviner, 'linear', 'extrap');
                        fprintf('\b (%g) \n', xx)
                    end
                else
                    new_value = interp1(T_model_interp_arr, value_interp_arr, T_diviner, 'linear', 'extrap');
                end
                if new_value == value_interp_arr(end)
                    if T_diviner < T_model
                        new_value = new_value*1.01; % ensure doesn't get stuck
                    elseif T_diviner > T_model
                        new_value = new_value*0.99;
                    end
                    if debug
                        fprintf('stuck value \t%i\t%i\t%g > %g\n', lat_idx, long_idx, value_interp_arr(2), new_value)
                    end
                end
                
                if new_value < value_limits(1)
                    new_value = value_limits(1);
                elseif new_value > value_limits(2)
                    new_value = value_limits(2);
                    if debug
                        fprintf('limit top\t%i\t%i\t%g > %g\n', lat_idx, long_idx, value_interp_arr(2), new_value)
                    end
                end
                
                if debug
                    if new_value > 0.5
                        fprintf('BIG: %i %i %g %g %g\n', lat_idx, long_idx, new_value, T_model, T_diviner)
                    end
                end
                value_3dmat(lat_idx, long_idx, loop_idx+1) = new_value;
            end
        end
    end
    
    if debug
        figure(1)
        plot_3d_surface(ls_data.ew_matrix, ls_data.ns_matrix, ls_data.elevation_matrix, value_3dmat(:,:,end), sprintf('Current values (after run %i/%i)', loop_idx, break_idx))
        figure(2)
        plot_3d_surface(ls_data.ew_matrix, ls_data.ns_matrix, ls_data.elevation_matrix, value_3dmat(:,:,end-1), 'Previous values')
        figure(3)
        plot_3d_surface(ls_data.ew_matrix, ls_data.ns_matrix, ls_data.elevation_matrix, value_3dmat(:,:,end) - value_3dmat(:,:,end-1), 'Value change')
        
        M = value_3dmat(:,:,end) - value_3dmat(:,:,end-1);
        if min(M(:)) - max(M(:)) == 0
            caxis(min(M(:) + 1e-3*[-1 1]))
        else
            caxis([min(M(:)), max(M(:))])
        end
        step_colormap(spectral, 1, 0, true);
        figure(4)
        plot_3d_surface(ls_data.ew_matrix, ls_data.ns_matrix, ls_data.elevation_matrix, T_model_3dmat(:,:,end), 'Current T')
        figure(5)
        plot_3d_surface(ls_data.ew_matrix, ls_data.ns_matrix, ls_data.elevation_matrix, T_model_3dmat(:,:,end) - T_diviner_matrix, 'Current T - Diviner T')
        divergent_colormap(10)
        drawnow
        fprintf('MAX: %g\n', max(max(value_3dmat(:,:,end))))
    end
end

data_filename = create_static_path('outputs/ls_parameters');
data_filename = sprintf('%s/%s_%s_parameter_fit', data_filename, ls_name, variable_name);
data = struct;
data.variable = variable_name;
data.value_matrix = value_3dmat(:,:,end);
data.ls_data = ls_data;
data.value_limits = value_limits;
data.break_idx = break_idx;
data.num_iterations = loop_idx;
data.fixed_data_paths = fixed_data_paths;
data.parameters = parameters;
data.T_mode = T_mode;
data.t_wait = t_wait;
data.filename = data_filename;
data.description = 'Matrix of parameter value for landing site, found by repeated interpolation comparing modelled and measured Diviner temperatures';
data.created = datetime;
save(data_filename, 'data')
fprintf('Parameters saved\n')

rehash
close all
clf
plot_ls_parameters(ls_name, variable_name, true);

if run_model
    rehash
    simulate_ls_temperatures(ls_name, 'fit');
end
end