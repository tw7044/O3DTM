function [output, theta_arr] = hayne_1d_model(lat, t_wait, T_mode, aspect, slope, custom_parameters, theta_arr, long, dtm_arr, sub_sun_long_arr, decl_arr)
% HAYNE_1D_MODEL is a 1D temperature model for the moon
%  Arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  lat = latitude of point in degrees
%  t_wait = controls how long the model runs for:
%  > number => number of days to run simulation for (to stabilise T) with
%    no seasonal variation
%  > 'seasons' => run with declination data to simulate lunar seasons
%  > datetime => run with seasons
%  T_mode = controls output of model:
%  > 'max' => output = maximum lunar surface temperature
%  > 'min' => output = minimum lunar surface temperature
%  > 'max time'/'min time' => output is struct containing the max/min
%    temperature and the time it occurs (fractions of lunar period since
%    midday)
%  > 'times' => output = struct containing arrays of times and their
%    corresponsing surface temperatures
%  > number (e.g. 0.5)  => output = array of layer temperatures at their
%    respective depths at this fraction of a lunar day after midday (e.g.
%    0.5 gives array of layer temperatures at midnight)
%  Optional Arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  aspect = slope direction in degrees (angle between horizontal projection
%  of surface normal and local meridian)
%  slope = slope magnitude in degrees (angle between surface normal and
%  vertical)
%  custom_parameters = struct of custom parameters for code
%  long = longitude for use with seasonal calculations (necessary for
%  seasonal calculations)
%  decl_arr = precomputed declination values (optional for seasonal
%  calculations - will speed up computation)
%  decl_dt_arr = precomputed declination times (optional for seasonal
%  calculations - will speed up computation)
%  theta_arr = precomputed solar incidence values (optional - will speed up
%  computation for seasonal calculations)
%
% This function is only dependent on define_parameters.m
% This model is based on the one described in Hayne et al. 2017
% By Oliver King, 2017

%% Process code arguments & load parameters
narginchk(3,inf)
if nargin < 4
    slope = 0;
    aspect = 0;
end
if nargin < 6
    custom_parameters = struct;
end
if nargin < 7 || numel(theta_arr) == 0
    compute_theta_arr = true;
else
    compute_theta_arr = false;
end
if nargin < 8
    long = 0;
else
    long = mod(long, 360);
end
if (ischar(t_wait) && strcmp(t_wait, 'seasons')) || isdatetime(t_wait)
    use_seasons = true;
else
    use_seasons = false;
end
if nargin < 9 && compute_theta_arr && use_seasons
    [dtm_arr, sub_sun_long_arr, decl_arr] = read_horizons_data();
end

parameters = define_parameters(custom_parameters);
% load parameters from define_parameters and convert to individual
% variables to increase code legibility and performance (e.g. using
% parameters.P is slower than using P)
P = parameters.P;
rho_s = parameters.rho_s;
rho_d = parameters.rho_d;
H = parameters.H;
stefans_constant = parameters.stefans_constant;
emissivity = parameters.emissivity;
S = parameters.S;
R_AU = parameters.R_AU;
Q = parameters.Q;
A0 = parameters.A0;
a = parameters.a;
b = parameters.b;
Chi = parameters.Chi;
Ks = parameters.Ks;
Kd = parameters.Kd;
c0 = parameters.c0;
c1 = parameters.c1;
c2 = parameters.c2;
c3 = parameters.c3;
c4 = parameters.c4;
m = parameters.m;
n = parameters.n;
num_skin_depths = parameters.num_skin_depths;
decl = parameters.decl;
T_min = parameters.T_min;
T_max = parameters.T_max;
T_bottom_limit = parameters.T_bottom_limit;
surface_bc_test_difference = parameters.surface_bc_test_difference;
surface_bc_break_counter = parameters.surface_bc_break_counter;
live_graph = parameters.live_graph;
live_graph_plot_interval = parameters.live_graph_plot_interval;
show_results_table = parameters.show_results_table;
return_z_grid = parameters.return_z_grid;
return_T_arr = parameters.return_T_arr;
dt = parameters.dt;
diviner_start_dtm = parameters.diviner_start_dtm;
diviner_end_dtm = parameters.diviner_end_dtm;
phase_start_dtm = parameters.phase_start_dtm;
phase_end_dtm = parameters.phase_end_dtm;
initial_depth_t_wait = parameters.initial_depth_t_wait;
grow_depth_t_wait = parameters.grow_depth_t_wait;
initial_num_skin_depths = parameters.initial_num_skin_depths;


% pre-calculate constant parameters
S_OVER_R_AU2 = S/R_AU^2;
emissivity_TIMES_stefans_constant = emissivity*stefans_constant;
Chi_OVER_350_POWER_3 = Chi*(1/350)^3;
a_OVER_pi_OVER_4_POWER_3 = a*(1/(pi/4))^3;
b_OVER_pi_OVER_2_POWER_8 = b*(1/(pi/2))^8;

%% Prepare declination data
% If using seasons, process declination data for use in simulation
use_ltim_end = false;
if use_seasons
    long_interp_arr = sub_sun_long_arr((dtm_arr > phase_start_dtm) & (dtm_arr < phase_end_dtm));
    dtm_interp_arr = dtm_arr((dtm_arr > phase_start_dtm) & (dtm_arr < phase_end_dtm));
    start_dtm = interp1(long_interp_arr(2:end-1), dtm_interp_arr(2:end-1), long); % start simulation at local midday
    if isnat(start_dtm)
        start_dtm = dtm_interp_arr(1); % start at sim at start of dtm range
    end
    
    if isdatetime(t_wait) && isnumeric(T_mode)
        % run to date then run to ltim
        measure_start_dtm = t_wait;
        end_dtm = t_wait + seconds(P*1.5);
        t_wait = seconds(t_wait - start_dtm)/P;
        use_ltim_end = true;
    elseif isdatetime(T_mode)
        % run to date
        end_dtm = T_mode;
        t_wait = 0;
    else
        % run over whole period (max/min)
        end_dtm = diviner_end_dtm;
        t_wait = seconds(diviner_start_dtm - start_dtm)/P;
    end
    
    t_limit = seconds(end_dtm - start_dtm);
    t_arr = start_dtm:seconds(dt):end_dtm;
    decl_arr = interp1(dtm_arr, decl_arr, t_arr);
    sub_sun_long_arr = interp1(dtm_arr, sub_sun_long_arr, t_arr);
    h_arr = long - sub_sun_long_arr;
    h_arr = mod(h_arr + 180, 360) - 180;
    
    if use_ltim_end
       % get appropriate end for ltim data
       h_end = T_mode*360;
       h_end = mod(h_end + 180, 360) - 180;
       idx = round(interp1(t_arr, 1:numel(t_arr), measure_start_dtm));
       h_end_big = h_end > max(h_arr);
       h_end_small = h_end < min(h_arr);
       for idx = idx:numel(h_arr)
           if (h_arr(idx-1) < h_end) && (h_arr(idx) >= h_end)
               % on rising slope
               break
           elseif (h_end_big) && (h_arr(idx-1) > h_arr(idx))
               % h_end = 180
               break
           elseif (h_end_small) && (h_arr(idx-1) > h_arr(idx))
               % h_end = -180
               break
           end
       end
       t_arr = t_arr(1:idx);
       decl_arr = decl_arr(1:idx);
       h_arr = h_arr(1:idx);
       t_limit = seconds(t_arr(end) - t_arr(1));
    end
else
    if ischar(T_mode)
        t_limit = P*(t_wait + 1) - dt;
    else
        t_limit = P*(t_wait + T_mode);
    end
end

%% Generate solar incidence angles
% Generate solar incidence angles using formula from Braun and Mitchell [1983]
if lat < 0
    aspect = -aspect;
    slope = -slope;
end
t_steps = P/dt;
dh = 360/t_steps;
if use_seasons
    t_steps = numel(decl_arr);
end
if compute_theta_arr
    theta_arr = NaN(1,t_steps);
    for t_idx = 1:t_steps
        if use_seasons
            decl = decl_arr(t_idx);
            h = h_arr(t_idx);
        else
            h = dh*t_idx - dh;
        end
        h = mod(h + 180, 360) - 180;
        if abs(h) < acosd(cotd(lat)*tand(decl))
            sigma_ew = 1;
        else
            sigma_ew = -1;
        end
        if lat*(lat-decl) >= 0
            sigma_ns = 1;
        else
            sigma_ns = -1;
        end
        if h >= 0
            sigma_w = 1;
        else
            sigma_w = -1;
        end
        theta_z = acos(sind(decl)*sind(lat) + cosd(decl)*cosd(lat)*cosd(h));
        if theta_z ~= 0
            gamma_so = real(asind(sind(h)*cosd(decl)/sin(theta_z)));
            % gamma_so may end up complex if it is ~90deg as rounding errors in
            % the asind(...) function can make its argument slightly larger
            % than 1, making gamma_so slightly complex. Therefore, only use the
            % real part of gamma_so to avoid the tiny complex part affecting
            % later calculations.
        else
            gamma_so = 0; % Avoid division by 0 error
        end
        gamma_s = sigma_ew*sigma_ns*gamma_so + ((1-sigma_ew*sigma_ns)/2)*sigma_w*180;
        theta = acos(cos(theta_z)*cosd(slope) + sin(theta_z)*sind(slope)*cosd(gamma_s - aspect));
        if theta < pi/2 && theta_z < pi/2
            % want sun to be visible from surface and to be above horizon
            theta_arr(t_idx) = theta;
        else
            theta_arr(t_idx) = NaN;
        end
    end
end


%% Generate grid
% Use smallest expected skin depth, zs_min to define thickness of layers
% and use largest expected skin depth, zs_max to define total number of
% layers to use
    function [K_value] = K_function(T, rho, par)        
        Kc_value = par.Kd - (par.Kd - par.Ks)*(par.rho_d - rho)/(par.rho_d - par.rho_s);
        K_value = Kc_value*(1 + par.Chi*(T/350)^3);
    end
    function [c_p_value] = c_p_function(T, par)
        c_p_value = par.c0 + par.c1*T + par.c2*T^2 + par.c3*T^3 + par.c4*T^4;
    end
zs_max = sqrt((K_function(T_max, rho_d, parameters)/(rho_s*c_p_function(T_max, parameters)))*P/pi);
zs_min = sqrt((K_function(T_min, rho_s, parameters)/(rho_s*c_p_function(T_min, parameters)))*P/pi);
layer_idx = 0;
dz = zs_min/m;
z = 0;
while z < zs_max*(num_skin_depths + 1)
    layer_idx = layer_idx + 1;
    z_arr(layer_idx) = z;
    dz_arr(layer_idx) = dz;
    if H == 0 && z == 0
        rho_arr(layer_idx) = rho_s;
    else
        rho_arr(layer_idx) = rho_d - (rho_d - rho_s)*exp(-z/H);
    end
    dz = dz*(1+1/n);
    z = z + dz;
end
num_layers = layer_idx;

if return_z_grid
    output = z_arr;
    return
end

initial_depth = zs_max*initial_num_skin_depths;
bottom_layer_idx = round(interp1(z_arr, 1:numel(z_arr), initial_depth));
depth_update_wait_t = initial_depth_t_wait;
depth_update_t_interval = grow_depth_t_wait/(num_layers - bottom_layer_idx - 1);


%% Initialise temperature profile
% Assume surface in eqilibrium with F_solar and black body emission, and
% deep temperature drops off to T0/sqrt(2)
% assumes t=0 start => hour angle = 0 at start etc.
theta = theta_arr(1);
albedo = A0 + a_OVER_pi_OVER_4_POWER_3*theta*theta*theta + b_OVER_pi_OVER_2_POWER_8*theta*theta*theta*theta*theta*theta*theta*theta;
T0 = (((1-albedo)/(emissivity*stefans_constant))*S/R_AU^2*cos(theta))^0.25;
if theta > pi/2 || T0 < T_bottom_limit || ~isreal(T0) || isnan(T0)
    % ensure that initial value doesn't break c_p etc.
    % try falling back to smooth moon, decl=0 temperature, and if not, just
    % set T0 to T_bottom_limit
    theta = abs(acos(cosd(lat)*cosd(0)*cosd(0) + sind(lat)*sind(0)));
    albedo = A0 + a_OVER_pi_OVER_4_POWER_3*theta*theta*theta + b_OVER_pi_OVER_2_POWER_8*theta*theta*theta*theta*theta*theta*theta*theta;
    T0 = (((1-albedo)/(emissivity*stefans_constant))*S/R_AU^2*(cos(theta)))^0.25;
    if T0 < T_bottom_limit || ~isreal(T0) || isnan(T0)
        T0 = T_bottom_limit;
    end
end
TN = T0/sqrt(2);
if TN < T_bottom_limit
    % set roughly constant initial temperature profile for very
    % cold regions
    TN = T_bottom_limit;
end
T_arr = NaN(1,num_layers);
for layer_idx = 1:num_layers
    z = z_arr(layer_idx);
    if z == 0 && H == 0
        T_arr(layer_idx) = T0;
    else
        T_arr(layer_idx) = TN - (TN - T0)*exp(-z/H);
    end
end

T_new_arr = NaN(1,num_layers);

%% Initialise other variables
if live_graph
    clf
    live_graph_idx = 1;
    T_layers_matrix = [];
    plotting_t_arr = [];
    graph_colors = colormap(viridis(num_layers));
end
if ischar(T_mode)
    Tsurface_arr = NaN(1,P/dt);
    time_arr = NaN(1,P/dt);
    Tsurface_idx = 1;
end

%% Prepare for loop
% Prepare for loop by pre-calculating constant values which do not need to
% be calculated every iteration (therefore increasing performance)
Kc_arr = NaN(1,num_layers);
for layer_idx = 1:num_layers
    Kc_arr(layer_idx) = Kd - (Kd - Ks)*(rho_d - rho_arr(layer_idx))/(rho_d - rho_s);
end

p_arr = NaN(1,num_layers);
q_arr = NaN(1,num_layers);
for layer_idx = 2:num_layers-1
    d3z = dz_arr(layer_idx)*dz_arr(layer_idx-1)*(dz_arr(layer_idx) + dz_arr(layer_idx-1));
    p = 2*dz_arr(layer_idx)/d3z;
    q = 2*dz_arr(layer_idx-1)/d3z;
    p_arr(layer_idx) = p;
    q_arr(layer_idx) = q;
end
B_surface = Kc_arr(1)*Chi_OVER_350_POWER_3;

%% Perform temperature simulation
% Loop through time steps performing temperature simulation
% Function calls take a relatively large amount of time, so functions e.g.
% K() are performed inline within this loop and exponents (e.g. T^n) are
% calculated by multiplication (T*T*T...) as this is faster in MATLAB
t_arr = 0:dt:t_limit;
for t_idx = 1:numel(t_arr)
    t = t_arr(t_idx);
    
    %% Increase simulation depth
    if bottom_layer_idx < num_layers && t > depth_update_wait_t
        depth_update_wait_t = depth_update_wait_t + depth_update_t_interval;
        % extend simulation by 1 layer, using bottom BC to define new
        % temperature
        bottom_layer_idx = bottom_layer_idx + 1;
        T_above = T_arr(bottom_layer_idx-1);
        T_arr(bottom_layer_idx) = T_above + dz_arr(bottom_layer_idx-1)*Q/(Kc_arr(bottom_layer_idx-1)*(1 + Chi_OVER_350_POWER_3*T_above*T_above*T_above));
    end
    
    %% Live graph plotting
    % Plot live graph if enabled (generally for debugging purposes as this
    % will decrease execution speed my many orders of magnitude)
    if live_graph
        cla
        T_layers_matrix(live_graph_idx,:) = T_arr;
        plotting_t_arr(end + 1) = t/P;
        if mod(live_graph_idx, live_graph_plot_interval) == 0 || t == 0 || t + dt > t_limit
            hold on
            for i = 1:num_layers
                plot(plotting_t_arr, T_layers_matrix(:,i), '.', 'Color', graph_colors(i,:), 'MarkerSize', 3);
            end
            if t == 0
                ylabel('Temperature (K)')
                xlabel('Lunar days')
                c = colorbar;
                c.Direction = 'reverse';
                c.Label.String = 'Layer depth (m)';
%                 c.YTick = z_arr;

                caxis([z_arr(1), z_arr(end)]);
                colormap(viridis)
            end
            drawnow
        end
        live_graph_idx = live_graph_idx + 1;
    end
    
    %% Surface BC
    break_iter = 0;
    difference = surface_bc_test_difference + 1;
    % Calculate incoming solar energy
    if use_seasons
        theta = theta_arr(t_idx);
    else
        h = 360*mod(t,P)/P; % hour angle
        theta = theta_arr(floor(h/dh + 1));
    end
    if isnan(theta)
        Qs = 0;
    else
        albedo = A0 + a_OVER_pi_OVER_4_POWER_3*theta*theta*theta + b_OVER_pi_OVER_2_POWER_8*theta*theta*theta*theta*theta*theta*theta*theta;
        Qs = (1 - albedo)*S_OVER_R_AU2*cos(theta);
    end
    T = T_arr(1);
    % Root finding
    while difference > surface_bc_test_difference
        % MY VERSION (from Hayne 2017)
        K_value = Kc_arr(1)*(1 + Chi_OVER_350_POWER_3*T*T*T);
        func = emissivity_TIMES_stefans_constant*T*T*T*T - Qs - K_value*((-3*T + 4*T_arr(2) - T_arr(3))/(2*dz_arr(1)));
        diff = 4*emissivity_TIMES_stefans_constant*T*T*T - 3*B_surface*T*T*((4*T_arr(2) - 3*T - T_arr(3))/(2*dz)) + (3/(2*dz_arr(1)))*K_value;
        
        % OLD VERSION:
%         func = emissivity_TIMES_stefans_constant*T*T*T*T - Qs - Kc_arr(1)*(4*T_arr(2) - 3*T - T_arr(3))/(2*dz_arr(1)) - B_surface*((T-T_arr(2))/2)^3;
%         diff = 4*emissivity_TIMES_stefans_constant*T*T*T + 3*Kc_arr(1)/(2*dz_arr(1)) - 3/2*B_surface*((T+T_arr(2))/2)^2*(T_arr(2)-T)/dz_arr(1);
        
        T_new = T - func/diff;
        if T_new <= 0
            fprintf('\nNegative T in surface BC, returning NaN')
            output = NaN;
            return
        end
        difference = abs(T_new - T);
        T = T_new;
        break_iter = break_iter + 1;
        if break_iter > surface_bc_break_counter
            fprintf('\nSurface BC iteration limit, returning NaN')
            output = NaN;
            return
        end
    end
    T_new_arr(1) = T_new;
    
    %% Bottom BC
    T_above = T_arr(bottom_layer_idx-1);
    T_new_arr(bottom_layer_idx) = T_above + dz_arr(bottom_layer_idx-1)*Q/(Kc_arr(bottom_layer_idx-1)*(1 + Chi_OVER_350_POWER_3*T_above*T_above*T_above));

    %% Middle layers
    for layer_idx = 2:bottom_layer_idx-1
        T = T_arr(layer_idx);
        T_above = T_arr(layer_idx-1);
        alpha_var = p_arr(layer_idx)*Kc_arr(layer_idx-1)*(1 + Chi_OVER_350_POWER_3*T_above*T_above*T_above);
        beta_var = q_arr(layer_idx)*Kc_arr(layer_idx)*(1 + Chi_OVER_350_POWER_3*T*T*T);
        c_p_value = c0 + c1*T + c2*T*T + c3*T*T*T + c4*T*T*T*T;
        T_new_arr(layer_idx) = T + (dt/(rho_arr(layer_idx)*c_p_value))*(alpha_var*T_above - (alpha_var+beta_var)*T + beta_var*T_arr(layer_idx+1));
    end
    
    %% Finish loop
    T_arr = T_new_arr;
    if t >= P*t_wait && ischar(T_mode)
        Tsurface_arr(Tsurface_idx) = T_arr(1);
        time_arr(Tsurface_idx) = t/P  - t_wait;
        Tsurface_idx = Tsurface_idx + 1;
    end
end
%% Prepare results & output
if show_results_table
    final_values_table = table(rot90(flip(z_arr)), rot90(flip(T_arr)));
    final_values_table.Properties.VariableNames = {'z','T_final'};
    disp(final_values_table)
end
if ischar(T_mode)
    if strcmp(T_mode, 'min')
        output = nanmin(Tsurface_arr);
    elseif strcmp(T_mode, 'max')
        output = nanmax(Tsurface_arr);
    elseif strcmp(T_mode, 'min time')
        [T, idx] = nanmin(Tsurface_arr);
        output = struct('T',T,'time', time_arr(idx));
    elseif strcmp(T_mode, 'max time')
        [T, idx] = nanmax(Tsurface_arr);
        output = struct('T',T,'time', time_arr(idx));
    else
        output = struct('T_arr', Tsurface_arr, 'time_arr', time_arr);
    end
    if use_seasons
        dtm_arr = start_dtm + seconds(P*(t_wait + time_arr));
        output.dtm_arr = dtm_arr;
    end
else
    if return_T_arr
        output = T_arr;
    else
        output = T_arr(1);
    end
end
end