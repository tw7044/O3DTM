% TEST_THETA_SUN debug routine to test solar angle ray tracing calculations

% close(figure(2))
ref_sphere = referenceSphere('moon');
r_moon = ref_sphere.Radius;

data = load_crater_environment('ls1_4ppd_wide_compressed');
lat_arr = data.lat_arr;% + 160;
long_arr = data.long_arr;
elevation_matrix = data.elevation_matrix-r_moon;

lat1_idx = 1;6;
long1_idx = 1;25;
lat2_idx = 3;
long2_idx = 10;

pause_seconds = 1;
plot_interval = 50;

lat1 = lat_arr(lat1_idx);
lat2 = lat_arr(lat2_idx);
long1 = long_arr(long1_idx);
long2 = long_arr(long2_idx);

figure(1)
plot_3d_surface(data.ew_matrix, data.ns_matrix, data.elevation_matrix)
c = colorbar;
delete(c);
hold on
scatter3(data.ew_matrix(lat1_idx,long1_idx), data.ns_matrix(lat1_idx,long1_idx), data.elevation_matrix(lat1_idx,long1_idx), 100, 'm.')
hh = plot3(NaN, NaN, NaN, 'b');
% scatter3(data.ew_matrix(lat2_idx,long2_idx), data.ns_matrix(lat2_idx,long2_idx), data.elevation_matrix(lat2_idx,long2_idx), 100, 'r.')

% ref = georefcells([min(lat_arr), max(lat_arr)], [min(long_arr), max(long_arr)], size(elevation_matrix), 'ColumnsStartFrom', 'north', 'RowsStartFrom', 'west');
% los2(elevation_matrix, ref, lat1, long1, lat2, long2, 0, 0, 'AGL', 'AGL', r_moon)


figure(2)
clf
figure(3)
clf
hold on
%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% 

parameters = define_parameters;
P = parameters.P;
const_decl = parameters.decl;
dt = parameters.dt;
phase_start_dtm = parameters.phase_start_dtm;
phase_end_dtm = parameters.phase_end_dtm;

long = mean(long_arr);
[raw_dtm_arr, sub_sun_long_arr, decl_arr] = read_horizons_data();
long_interp_arr = sub_sun_long_arr((raw_dtm_arr > phase_start_dtm) & (raw_dtm_arr < phase_end_dtm));
dtm_interp_arr = raw_dtm_arr((raw_dtm_arr > phase_start_dtm) & (raw_dtm_arr < phase_end_dtm));
% ensure both have consistent ranges
long = mod(long + 180, 360) - 180;
long_interp_arr = mod(long_interp_arr + 180, 360) - 180;
[long_interp_arr, unique_idx_arr] = unique(long_interp_arr);
dtm_interp_arr = dtm_interp_arr(unique_idx_arr);
start_dtm = interp1(long_interp_arr(2:end-1), dtm_interp_arr(2:end-1), long); % start simulation at local midday
end_dtm = raw_dtm_arr(end);
% start_dtm = datetime(2014,11,4,13,0,30);
dtm_arr = start_dtm:seconds(dt):end_dtm;
decl_arr = interp1(raw_dtm_arr, decl_arr, dtm_arr);
sub_sun_long_arr = interp1(raw_dtm_arr, sub_sun_long_arr, dtm_arr);
h_arr = long - sub_sun_long_arr;
h_arr = mod(h_arr + 180, 360) - 180;



extended_elevation_matrix = zeros(size(elevation_matrix) + [2,2]);
extended_elevation_matrix(2:end-1,2:end-1) = elevation_matrix;
extended_elevation_matrix(2:end-1,1) = elevation_matrix(:,1);
extended_elevation_matrix(2:end-1,end) = elevation_matrix(:,end);
extended_elevation_matrix(1,2:end-1) = elevation_matrix(1,:);
extended_elevation_matrix(end,2:end-1) = elevation_matrix(end,:);
extended_elevation_matrix(1,1) = elevation_matrix(1,1);
extended_elevation_matrix(1,end) = elevation_matrix(1,end);
extended_elevation_matrix(end,1) = elevation_matrix(end,1);
extended_elevation_matrix(end,end) = elevation_matrix(end,end);

d_lat = abs(lat_arr(1) - lat_arr(2));
d_long = abs(long_arr(1) - long_arr(2));
extended_raster_ref = georefcells([min(lat_arr)-d_lat, max(lat_arr)+d_lat], [min(long_arr)-d_long, max(long_arr)+d_long], size(extended_elevation_matrix), 'ColumnsStartFrom', 'north', 'RowsStartFrom', 'west');



[aspect_matrix, slope_matrix] = convert_height_to_slope(elevation_matrix, lat_arr, long_arr);
[~,~,max_distance] = geodetic2aer(...
    lat_arr(1), long_arr(1), min(reshape(elevation_matrix,[],1)),...
    lat_arr(end), long_arr(end), max(reshape(elevation_matrix,[],1)),...
    ref_sphere);

lat_idx = lat1_idx;
long_idx = long1_idx;
long = long_arr(long_idx);
lat = lat_arr(lat_idx);



%% Get information about current location
slope1 = slope_matrix(lat_idx, long_idx);
aspect1 = aspect_matrix(lat_idx, long_idx);
height_px = elevation_matrix(lat_idx, long_idx);
if lat >= 0
    aspect = aspect1;
    slope = slope1;
else
    aspect = -aspect1;
    slope = -slope1;
end

h_arr_local = h_arr + long - mean(long_arr);
decl_arr_local = decl_arr;


%% Generate visibility matrix
% vis_4dmat_local(1, long_idx, :, :) = viewshed(elevation_matrix, raster_ref, lat, long, 0, 0, 'AGL', 'AGL', r_moon);

%% Generate solar incidence angles & solar fluxes
% Generate solar incidence angles using formula from Braun and Mitchell [1983]
t_steps = numel(decl_arr);
fprintf('azi.\telev.\ttta_z \ttta_s \tvis\n')

for t_idx = 1:t_steps
    decl = decl_arr_local(t_idx);
    h = h_arr_local(t_idx);
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
    
    if true || (theta < pi/2 && theta_z < pi/2)
        % want sun to be visible from surface and to be above horizon
        if theta_z == 0
            sun_vis = true;
        else
            elev_sun = 90-rad2deg(theta_z);
            if lat >= decl
                az_sun = 180 + gamma_s;
            else
                az_sun = -gamma_s;
            end
            if mod(t_idx-1,plot_interval) == 0
                slant_range = 1.1*max_distance; % ensure outside of grid
                [lat_sun,long_sun,height_sun] = aer2geodetic(az_sun,elev_sun,slant_range,lat,long,height_px,ref_sphere);
                [sun_vis,visprofile,dist,H,lattrk,lontrk] = los2(extended_elevation_matrix, extended_raster_ref, lat, long, lat_sun, long_sun, 0, height_sun, 'AGL', 'MSL', r_moon);
                fprintf('%3.1f \t%3.1f \t%3.1f \t%3.1f \t%g\n',az_sun, elev_sun, rad2deg(theta_z), rad2deg(theta), sun_vis)
                
                Z = extended_elevation_matrix;
                R = extended_raster_ref;
                lat1 = lat;
                lon1 = long;
                oalt = 0;
                lat2 = lat_sun;
                lon2 = long_sun;
                talt = height_sun;
                observerAltitudeIsAGL = true;
                targetAltitudeIsAGL = false;
                actualradius = r_moon;
                F = griddedInterpolant(Z);
                apparentradius = actualradius;
                
                
                [visprofile, dist, h, lattrk, lontrk, x1, z1, x2, z2] = calculateLOS(F, R, ...
                    lat1, lon1, lat2, lon2, oalt, talt, observerAltitudeIsAGL, ...
                    targetAltitudeIsAGL, actualradius, apparentradius);
                
                figure(2)
                cla
                plotProfile(x1, z1, x2, z2, visprofile)
                title('Line of sight to sun')
                
                
                figure(1)
                [xEast,yNorth,zUp] = aer2enu(az_sun,elev_sun,slant_range);
                x = [0 xEast] + data.ew_matrix(lat1_idx,long1_idx);
                y = [0 yNorth] + data.ns_matrix(lat1_idx,long1_idx);
                z = [0 zUp] +data.elevation_matrix(lat1_idx,long1_idx);
                delete(hh)
                if sun_vis
                    c_str = 'b';
                else
                    c_str = 'r';
                end
                hh = plot3(x, y, z, c_str);
                title('Direction to sun')
                
                
                figure(3)
                scatter(t_idx, 90-rad2deg(theta), 100, sprintf('%s.',c_str))
                xlabel('Time step')
                ylabel('\theta_{sun} (degrees above parallel to surface)')
                title('Theta sun variation')
                grid on
                drawnow
                pause(pause_seconds)
            end
        end
    end
end
    
    
function plotProfile(x, z, x2, z2, vis)
% Plot the terrain profile plus observer position, visible and obscured
% points, and line of sight to last point in profile, in a new figure with
% the tag 'los2'.

vis = reshape(vis, size(x));
% 
% figure('Tag','los2','NextPlot','add')
% ax = axes('NextPlot','add');
ax = gca;
hold on
h(1) = plot(ax,x,z,'k');

if vis(end)
    h(5) = plot(ax,[0;x2(end)],[0;z2(end)],'b');
else
    h(5) = plot(ax,[0;x2(end)],[0;z2(end)],'r');
end

if any(vis)
    h(2) = plot(ax,x2(vis),z2(vis),'b+');
end

if any(~vis)
    h(3) = plot(ax,x2(~vis),z2(~vis),'ro');
end

h(4) = plot(ax,0,0,'mx');
axis(ax,'equal')

labels = {'Terrain','Visible','Obscured','Observer','Line of Sight'};
indx = ishghandle(h,'line');
legend(ax, h(indx), labels(indx))
xlabel(ax,'Horizontal Distance from Observer')
ylabel(ax,'Vertical Distance from Observer')
end