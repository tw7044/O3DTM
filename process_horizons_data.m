function [dtm_arr, sub_sun_long_arr, decl_arr] = process_horizons_data(plot_graph)
% PROCESS_HORIZONS_DATA converts horizons ephemeris .txt outputs into
% MATLAB matrices
%
% [dtm_arr, sub_sun_long_arr, decl_arr] = PROCESS_HORIZONS_DATA loads .txt
% files and returns matlab matrices of data
%
% ... = PROCESS_HORIZONS_DATA(plot_graph) controls if graph should be
% plotted of solar values

raw_data_path_root = create_static_path('ephemerides/horizons_results_');
raw_data_path_end = '.txt';
start_year = 2000;
end_year = 2015;

dtm_arr = {};
sub_sun_long_arr = {};
decl_arr = {};

%% Read in data
for year = start_year:end_year
    raw_data_path = sprintf('%s%i%s', raw_data_path_root, year, raw_data_path_end);
    fid = fopen(raw_data_path);
    if fid == -1
        error('File "%s" not found', raw_data_path)
    end
    header_line = true;
    while ~feof(fid)
        file_line = fgetl(fid);
        if strcmp(file_line, '$$SOE')
            header_line = false;
        elseif header_line
            continue
        elseif strcmp(file_line, '$$EOE')
            break
        else
            dtm_arr{end+1} = file_line(2:18);
            sub_sun_long_arr{end+1} = file_line(32:40);
            decl_arr{end+1} = file_line(42:50);
        end
    end
    fclose(fid);
end
dtm_arr = datetime(dtm_arr);
sub_sun_long_arr = str2double(sub_sun_long_arr);
decl_arr = str2double(decl_arr);

%% Ensure sub_sun_long loops smoothly 0<->360
% Make long change as ~step function as it goes from 0 to 360 so that there
% aren't issues when the values are interpolated later
idx = 2;
cycle_idx = 1;
while true
    if sub_sun_long_arr(idx) > sub_sun_long_arr(idx - 1)
        dtm1 = interp1(...
            [sub_sun_long_arr(cycle_idx+1:idx), sub_sun_long_arr(idx+1)-360],...
            dtm_arr(cycle_idx+1:idx+1),...
            0);
        long1 = 0;
        decl1 = interp1(dtm_arr, decl_arr, dtm1);
        dtm2 = dtm1 + seconds(1e-3);
        long2 = 360;
        decl2= interp1(dtm_arr, decl_arr, dtm2);
        dtm_arr = [dtm_arr(1:idx-1), dtm1, dtm2, dtm_arr(idx:end)];
        sub_sun_long_arr = [sub_sun_long_arr(1:idx-1), long1, long2, sub_sun_long_arr(idx:end)];
        decl_arr = [decl_arr(1:idx-1), decl1, decl2, decl_arr(idx:end)];
        idx = idx + 2;
        cycle_idx = idx;
    end
    idx = idx + 1;
    if idx > numel(dtm_arr)
        break
    end
end

if nargin > 0 && plot_graph
    clf
    plot(dtm_arr, sub_sun_long_arr)
    hold on
    yyaxis right
    plot(dtm_arr, decl_arr)
    legend('Sub sun long', 'Declination')
    hold off
end
end