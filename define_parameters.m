function parameters = define_parameters(custom_parameters)
% DEFINE_PARAMETERS returns structure containing model parameters
%
% parameters = DEFINE_PARAMETERS returns the default model parameters
% 
% parameters = DEFINE_PARAMETERS(custom_parameters) returns model
% parameters with default values replaced by custom parameter values,
% defined by struct custom_parameters. E.g. custom_parameters =
% struct('A0', 0.1) returns struct of default parameters with A0 value
% replaced with 0.1

parameters = struct;

%% Define default parameters (from table A1 in Hayne et al. 2017)
% all constants in SI units
parameters.P = 2.55024e6; % lunar diurnal period, s
parameters.rho_s = 1100; % surface layer density, kg m^-3  
parameters.rho_d = 1800; % deep layer density, kg m^-3
parameters.H = 0.06;
parameters.stefans_constant = 5.6704e-8; % from data sheet
parameters.emissivity = 0.95; % IR emissivity
parameters.S = 1361; % solar constant, W m^-2
parameters.R_AU = 1; % astronomical unit **in AU**
parameters.Q = 0.018; % interior heat flow, W m^-2
parameters.A0 = 0.12; % bond albedo at normal solar incedence

% parameters.a = 0.06; % coefficinet for calculating albedo
% parameters.b = 0.25; % coefficient for calculating albedo
parameters.a = 0.0162; % coefficient for calculating albedo (new value from Apollo 16)
parameters.b = -0.03625; % coefficient for calculating albedo (new value from Apollo 16)

parameters.Chi = 2.7; % radiative condictivity parameter
parameters.Ks = 7.4e-4; % surface layer conductivity, W m^-1 K^-1
parameters.Kd = 3.4e-3; % deep layer conductivity, W m^-1 K^-1
parameters.c0 = -3.6125; % Coefficients for specific heat capacity fn, J kg^-1 K^-1
parameters.c1 = 2.7431;
parameters.c2 = 2.3616e-3;
parameters.c3 = -1.2340e-5;
parameters.c4 = 8.9093e-9;

%% Define dates (all UTC)
% parameters.phase_start_dtm = datetime(2004,1,7,8,45,0); % datetime for start of lunar phase cycle
% parameters.phase_end_dtm = datetime(2004,2,5,23,25,0); % datetime for end of lunar phase cycle
parameters.phase_start_dtm = datetime(2000,1,21,9,45,0); % datetime for start of lunar phase cycle
parameters.phase_end_dtm = datetime(2000,2,19,23,55,0); % datetime for end of lunar phase cycle
parameters.diviner_start_dtm = datetime(2009,7,5,16,0,0); % datetime for start of temperature data
parameters.diviner_end_dtm = datetime(2015,11,4,19,0,0); % datetime for end of temperature data


%% Define code controls
parameters.m = 10;
parameters.n = 5;
parameters.num_skin_depths = 10;
parameters.max_depth = 2.5; % depth to run 3d simulaton to (m)
parameters.initial_num_skin_depths = 2; % shallow depth for start of simulation run
parameters.initial_depth_t_wait = 1*(60*60*24*365); % run simulation with shallow depth at start for 1 year
parameters.grow_depth_t_wait = 3*(60*60*24*365); % run simulation with increasing depth for next 3 years
parameters.decl = 0; %deg
parameters.T_min = 93;
parameters.T_max = 388;
parameters.T_bottom_limit = 30; % approx. bottom limit of temperature from core heat flux to set floor for initialised temperature
parameters.surface_bc_test_difference = 1e-5; % accuracy to check convergence to for surface bc
parameters.surface_bc_break_counter = 1e3; % max. number of iterations for surface bc
parameters.live_graph = false; % adds massive overhead
parameters.live_graph_plot_interval = 10; % how many intervals between updating the live graph
parameters.show_results_table = false;
parameters.return_z_grid = false; % return z grid instead of data
parameters.return_T_arr = false; % return array of layer temperatures rather than just surface T
parameters.dt = parameters.P/(24*60*2); % time interval of 1/2 lunar minute
parameters.time_series_interval = 10; % number of dt time steps to save time series data


%% Update custom parameters
if nargin > 0
    parameters_to_update = fieldnames(custom_parameters);
    for i = 1:length(parameters_to_update)
        parameters.(parameters_to_update{i}) = custom_parameters.(parameters_to_update{i});
    end
end

%% Check parameters for issues

% Ensure Albedo is never greater than 1
if numel(parameters.A0) == 1 && parameters.A0 + parameters.a*2^3 + parameters.b > 1
        A0_new = 1 - parameters.a*2^3 - parameters.b;
        fprintf('Reducing A0 from %g to %g to ensure A <= 1 for all angles\n', parameters.A0, A0_new)
        parameters.A0 = A0_new;
elseif numel(parameters.A0) > 1
    for idx = 1:numel(parameters.A0)
        if parameters.A0(idx) + parameters.a*2^3 + parameters.b > 1
            A0_new = 1 - parameters.a*2^3 - parameters.b;
            fprintf('Reducing A0 from %g to %g to ensure A <= 1 for all angles\n', parameters.A0, A0_new)
            parameters.A0(idx) = A0_new;
        end
    end
end
    
if mod(parameters.P, parameters.dt) ~= 0
    dt_new = parameters.P/ceil(parameters.P/parameters.dt);
    fprintf('Changing dt from %g to %g to ensure integer time steps in P\n', parameters.dt, dt_new)
    parameters.dt = dt_new;
end

parameters.parameter_generation_datetime = datetime;
    
end