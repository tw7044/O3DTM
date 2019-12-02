% SAVE_HORIZONS_DATA saves processed JPL horizons data to a data file for
% use in other scripts 

fprintf('Processing horizons data...\n')
[dtm_arr, sub_sun_long_arr, decl_arr] = process_horizons_data;

fprintf('Saving data...\n')
horizons_data = struct;
horizons_data.dtm_arr = dtm_arr;
horizons_data.sub_sun_long_arr = sub_sun_long_arr;
horizons_data.decl_arr = decl_arr;
horizons_data.created = datetime;
horizons_data.description = 'Ephemeris data from JPL Horizons for lunar declination and the sub-sun longitude';

target_name = create_static_path('ephemerides/horizons_data.mat');
save(target_name, 'horizons_data');
fprintf('Done\n')