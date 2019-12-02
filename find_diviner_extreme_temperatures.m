function [Tmax_matrix, Tmin_matrix, Tmax_ltim_matrix, Tmin_ltim_matrix, Tmax_dtm_matrix, Tmin_dtm_matrix, lat_arr, long_arr] = ...
    find_diviner_extreme_temperatures(ppd, required_lat_arr, required_long_arr)
% FIND_DIVINER_EXTREME_TEMPERATURES finds extreme diviner temperatures for
% given resolution
%
% [Tmax_matrix, Tmin_matrix, Tmax_ltim_matrix, Tmin_ltim_matrix,
% Tmax_dtm_matrix, Tmin_dtm_matrix, lat_arr, long_arr] =
% FIND_DIVINER_EXTREME_TEMPERATURES(...) returns matrices of information
% about extreme temperatures
%
% ... = FIND_DIVINER_EXTREME_TEMPERATURES(ppd) finds diviner extreme
% temperatures for given resolution

% ... = FIND_DIVINER_EXTREME_TEMPERATURES(ppd, required_lat_arr,
% required_long_arr) finds diviner extreme temperatures for given
% latitude and longitude arrays

%% Set up
if nargin == 0
    ppd = 4;
end
switch ppd
    case 4
        tbol_file_location = create_static_path('4_diviner_tbol/');
        ltim_file_location = create_static_path('4_diviner_ltim/');
    case 16
        tbol_file_location = create_static_path('16_diviner_tbol/');
        ltim_file_location = create_static_path('16_diviner_ltim/');
    case 128
        tbol_file_location = create_static_path('128_diviner_tbol/');
        ltim_file_location = create_static_path('128_diviner_ltim/');
    otherwise
        error('Unexpected ppd')
end

tbol_jp2_names = dir([tbol_file_location,'*.jp2']);
ltim_jp2_names = dir([ltim_file_location,'*.jp2']);

if ppd == 128 || nargin == 3
    size_arr = [numel(required_lat_arr), numel(required_long_arr)];
    Tmax_matrix = ones(size_arr)*-inf;
    Tmin_matrix = ones(size_arr)*inf;
    Tmax_ltim_matrix = NaN(size_arr);
    Tmin_ltim_matrix = NaN(size_arr);
    Tmax_dtm_matrix = NaT(size_arr);
    Tmin_dtm_matrix = NaT(size_arr);
else
    Tmax_matrix = ones(180*ppd, 360*ppd)*-inf;
    Tmin_matrix = ones(180*ppd, 360*ppd)*inf;
    Tmax_ltim_matrix = NaN(180*ppd, 360*ppd);
    Tmin_ltim_matrix = NaN(180*ppd, 360*ppd);
    Tmax_dtm_matrix = NaT(180*ppd, 360*ppd);
    Tmin_dtm_matrix = NaT(180*ppd, 360*ppd);
end

%% Process files
fprintf('Finding extreme temperatures...\n')
progress_limit = 50;
fprintf('%s\n', repmat('_',1,progress_limit))
progress_idx = 0;
for file_idx = 1:length(tbol_jp2_names)
    %% Get filenames
    tbol_filename = tbol_jp2_names(file_idx).name;
    tbol_filename = tbol_filename(1:end-4);
    % find corresponding ltim filename
    ltim_filename = 'error';
    for ltim_file_idx = 1:length(ltim_jp2_names)
        if numel(strfind(ltim_jp2_names(ltim_file_idx).name, tbol_filename(19:27))) > 0
            ltim_filename = ltim_jp2_names(ltim_file_idx).name;
            ltim_filename = ltim_filename(1:end-4);
        end
    end
    if strcmp('error', ltim_filename)
        error('Cannot find ltim file for "%s"', tbol_filename)
    end
    
    %% Read ltim data
    fileID = fopen([ltim_file_location, ltim_filename, '.lbl'],'r');
    text_data = fscanf(fileID,'%s');
    
    pointer = strfind(text_data,'MAP_RESOLUTION');
    MAP_RESOLUTION = find_str(text_data,pointer);
    pointer = strfind(text_data,'CENTER_LONGITUDE');
    CENTER_LONGITUDE = find_str(text_data,pointer);
    pointer = strfind(text_data,'CENTER_LATITUDE');
    CENTER_LATITUDE = find_str(text_data,pointer);
    pointer = strfind(text_data,'LINE_PROJECTION_OFFSET');
    LINE_PROJECTION_OFFSET = find_str(text_data,pointer);
    pointer = strfind(text_data,'SAMPLE_PROJECTION_OFFSET');
    SAMPLE_PROJECTION_OFFSET = find_str(text_data,pointer);
    pointer = strfind(text_data,'MAXIMUM_LATITUDE');
    MAXIMUM_LATITUDE = find_str(text_data,pointer);
    pointer = strfind(text_data,'MINIMUM_LATITUDE');
    MINIMUM_LATITUDE = find_str(text_data,pointer);
    pointer = strfind(text_data,'WESTERNMOST_LONGITUDE');
    WESTERNMOST_LONGITUDE = find_str(text_data,pointer);
    pointer = strfind(text_data,'EASTERNMOST_LONGITUDE');
    EASTERNMOST_LONGITUDE = find_str(text_data,pointer);
    pointer = strfind(text_data,'A_AXIS_RADIUS');
    A_AXIS_RADIUS = find_str(text_data,pointer);
    pointer = strfind(text_data,'MAP_SCALE');
    MAP_SCALE = find_str(text_data,pointer);
    pointer = strfind(text_data,'LINE_LAST_PIXEL');
    LINE_LAST_PIXEL = find_str(text_data,pointer);
    pointer = strfind(text_data,'SAMPLE_LAST_PIXEL');
    SAMPLE_LAST_PIXEL = find_str(text_data,pointer);
    pointer = strfind(text_data,'OFFSET');
    OFFSET = find_str(text_data,pointer(2));
    pointer = strfind(text_data,'SCALING_FACTOR');
    SCALING_FACTOR = find_str(text_data,pointer(3));
    pointer = strfind(text_data,'MISSING_CONSTANT');
    MISSING_CONSTANT = find_str(text_data,pointer);
    pointer = strfind(text_data,'START_TIME');
    START_TIME = find_dt(text_data,pointer);
    pointer = strfind(text_data,'STOP_TIME');
    STOP_TIME = find_dt(text_data,pointer);
    pointer = strfind(text_data,'PRODUCT_VERSION_ID');
    PRODUCT_VERSION_ID = find_text(text_data,pointer);
    fclose(fileID);
    
    if file_idx > 1 && (ppd == 128 || nargin == 3)
        % save memory for 128ppd stuff
        ltim_matrix = imread([ltim_file_location,ltim_filename,'.jp2']);
        ltim_matrix = double(ltim_matrix(bool_matrix));
        ltim_matrix = reshape(ltim_matrix, size_arr);
    else
        ltim_matrix = double(imread([ltim_file_location,ltim_filename,'.jp2']));
    end
    ltim_matrix = (ltim_matrix*SCALING_FACTOR)+OFFSET;
    ltim_matrix(ltim_matrix == ((MISSING_CONSTANT*SCALING_FACTOR)+OFFSET)) = NaN;
    
    if file_idx == 1
        for LINE = 1:size(ltim_matrix,1)
            lat_arr(LINE,1) = CENTER_LATITUDE - (LINE - LINE_PROJECTION_OFFSET - 1) / MAP_RESOLUTION;
        end
        for SAMPLE = 1:size(ltim_matrix,2)
            long_arr(SAMPLE,1) = CENTER_LONGITUDE + (SAMPLE - SAMPLE_PROJECTION_OFFSET - 1) / MAP_RESOLUTION;
        end
        
        if ppd == 128 || nargin == 3
            lat_bool_arr = ismember(lat_arr, required_lat_arr);
            long_bool_arr = ismember(long_arr, required_long_arr);
            bool_matrix = reshape(lat_bool_arr, [], 1) & reshape(long_bool_arr, 1, []);
            lat_arr = lat_arr(lat_bool_arr);
            long_arr = long_arr(long_bool_arr);
            ltim_matrix = ltim_matrix(bool_matrix);
            ltim_matrix = reshape(ltim_matrix, size_arr);
        end
    end

    
    %% Generate datetimes
    dtm_matrix = repmat(START_TIME, size(ltim_matrix));

    %% Read tbol data
    fileID = fopen([tbol_file_location, tbol_filename,'.lbl'],'r');
    text_data = fscanf(fileID,'%s');
    
    pointer = strfind(text_data,'MAP_RESOLUTION');
    MAP_RESOLUTION = find_str(text_data,pointer);
    pointer = strfind(text_data,'CENTER_LONGITUDE');
    CENTER_LONGITUDE = find_str(text_data,pointer);
    pointer = strfind(text_data,'CENTER_LATITUDE');
    CENTER_LATITUDE = find_str(text_data,pointer);
    pointer = strfind(text_data,'LINE_PROJECTION_OFFSET');
    LINE_PROJECTION_OFFSET = find_str(text_data,pointer);
    pointer = strfind(text_data,'SAMPLE_PROJECTION_OFFSET');
    SAMPLE_PROJECTION_OFFSET = find_str(text_data,pointer);    
    pointer = strfind(text_data,'MAXIMUM_LATITUDE');
    MAXIMUM_LATITUDE = find_str(text_data,pointer);    
    pointer = strfind(text_data,'MINIMUM_LATITUDE');
    MINIMUM_LATITUDE = find_str(text_data,pointer);    
    pointer = strfind(text_data,'WESTERNMOST_LONGITUDE');
    WESTERNMOST_LONGITUDE = find_str(text_data,pointer);    
    pointer = strfind(text_data,'EASTERNMOST_LONGITUDE');
    EASTERNMOST_LONGITUDE = find_str(text_data,pointer);    
    pointer = strfind(text_data,'A_AXIS_RADIUS');
    A_AXIS_RADIUS = find_str(text_data,pointer);    
    pointer = strfind(text_data,'MAP_SCALE');
    MAP_SCALE = find_str(text_data,pointer);    
    pointer = strfind(text_data,'LINE_LAST_PIXEL');
    LINE_LAST_PIXEL = find_str(text_data,pointer);    
    pointer = strfind(text_data,'SAMPLE_LAST_PIXEL');
    SAMPLE_LAST_PIXEL = find_str(text_data,pointer);    
    pointer = strfind(text_data,'OFFSET');
    OFFSET = find_str(text_data,pointer(2));    
    pointer = strfind(text_data,'SCALING_FACTOR');
    SCALING_FACTOR = find_str(text_data,pointer(3));
    pointer = strfind(text_data,'MISSING_CONSTANT');
    MISSING_CONSTANT = find_str(text_data,pointer);
    
    fclose(fileID);
    if ppd == 128 || nargin == 3
        tbol_matrix = imread([tbol_file_location, tbol_filename,'.jp2']);
        tbol_matrix = double(tbol_matrix(bool_matrix));
        tbol_matrix = reshape(tbol_matrix, size_arr);
    else
        tbol_matrix = double(imread([tbol_file_location, tbol_filename,'.jp2']));
    end
    tbol_matrix = (tbol_matrix*SCALING_FACTOR)+OFFSET;
    tbol_matrix(tbol_matrix == ((MISSING_CONSTANT*SCALING_FACTOR)+OFFSET)) = NaN;
    tbol_matrix(tbol_matrix <= 0) = NaN; % remove obviously unphysical data points
    tbol_matrix(isnan(ltim_matrix)) = NaN;
    
    %% Process data
    Tmax_ltim_matrix(tbol_matrix > Tmax_matrix) = ltim_matrix(tbol_matrix > Tmax_matrix);
    Tmin_ltim_matrix(tbol_matrix < Tmin_matrix) = ltim_matrix(tbol_matrix < Tmin_matrix);
    
    Tmax_dtm_matrix(tbol_matrix > Tmax_matrix) = dtm_matrix(tbol_matrix > Tmax_matrix);
    Tmin_dtm_matrix(tbol_matrix < Tmin_matrix) = dtm_matrix(tbol_matrix < Tmin_matrix);

    Tmax_matrix(tbol_matrix > Tmax_matrix) = tbol_matrix(tbol_matrix > Tmax_matrix);
    Tmin_matrix(tbol_matrix < Tmin_matrix) = tbol_matrix(tbol_matrix < Tmin_matrix);

    if file_idx > progress_idx*(length(tbol_jp2_names)/progress_limit)
        fprintf('#')
        progress_idx = progress_idx + 1;
%         plot_map(lat_arr, long_arr, Tmax_ltim_matrix);
    end
end
fprintf('\n')
Tmax_matrix(isinf(Tmax_matrix)) = NaN;
Tmin_matrix(isinf(Tmin_matrix)) = NaN;
Tmax_ltim_matrix = mod(Tmax_ltim_matrix,24);
Tmin_ltim_matrix = mod(Tmin_ltim_matrix,24);

%% Define functions
    
    function back = find_str(text_data,pointer)
        for i=pointer:pointer+50
            if(text_data(i) == '=')
                pointer =i;
                break
            end
        end
        count = 1;
        for i=pointer:pointer+10
            if((double(text_data(i)) < 58 && double(text_data(i)) > 47) || double(text_data(i)) == 46 || double(text_data(i)) == 45)
                text_temp(count) = text_data(i);
                count = count +1;
            end
        end
        back = str2num(text_temp);
    end
    function back = find_text(text_data,pointer)
        for i=pointer:pointer+50
            if(text_data(i) == '"')
                pointer = i+1;
                break
            end
        end
        for i = pointer:pointer+100
            if(text_data(i) == '"')
                pointer_end = i-1;
                break
            end
        end
        back = text_data(pointer:pointer_end);
    end
    function back = find_dt(text_data,pointer)
        for i=pointer:pointer+50
            if(text_data(i) == '=')
                pointer = i;
                break
            end
        end
        text_data = text_data(pointer+1:pointer+23);
        back = datetime(text_data, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.S');
    end
end