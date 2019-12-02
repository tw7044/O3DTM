function output = save_ls_diviner_temperatures(ls_name)
% SAVE_LS_DIVINER_TEMPERATURES saves all diviner temperature measurements
% from a landing site into a data file
%
% output = SAVE_LS_DIVINER_TEMPERATURES(ls_name) returns copy of saved data
% file

%% Process inputs
if strcmp(ls_name(1:3), 'ls_')
    for ls_idx = flip(1:8)
        save_ls_diviner_temperatures(sprintf('ls%d%s', ls_idx, ls_name(3:end)));
    end
    return
end

ls_data = load_crater_environment(ls_name);
required_lat_arr = ls_data.lat_arr;
required_long_arr = ls_data.long_arr;
ppd = ls_data.ppd;

%% Uncompress if needed and prepare lat/long stuff
if isfield(ls_data, 'compression_ratio')
    compressed_environment = true;
    compressed_long_arr = required_long_arr;
    compression_ratio = ls_data.compression_ratio;
    compressed_ppd = ls_data.compressed_ppd;
    offset = 0.5*(1/compressed_ppd)*(compression_ratio-1)/(compression_ratio);
    offset_arr = linspace(-offset, offset, compression_ratio);
    required_long_arr = [];
    for compressed_long = compressed_long_arr
        required_long_arr(end+1:end+numel(offset_arr)) = compressed_long + offset_arr;
    end
else
    compressed_environment = false;
end
    
[long_matrix, lat_matrix] = meshgrid(required_long_arr, required_lat_arr);


%% Read in data file locations & set up
switch ppd
    case 4
        tbol_file_location = create_static_path('4_diviner_tbol/');
        ltim_file_location = create_static_path('4_diviner_ltim/');
        jd_file_location = create_static_path('4_diviner_jd/');
    case 16
        tbol_file_location = create_static_path('16_diviner_tbol/');
        ltim_file_location = create_static_path('16_diviner_ltim/');
        jd_file_location = create_static_path('16_diviner_jd/');

    case 128
        tbol_file_location = create_static_path('128_diviner_tbol/');
        ltim_file_location = create_static_path('128_diviner_ltim/');
        jd_file_location = create_static_path('128_diviner_jd/');
    otherwise
        error('Unexpected ppd')
end

tbol_jp2_names = dir([tbol_file_location,'*.jp2']);
ltim_jp2_names = dir([ltim_file_location,'*.jp2']);
jd_jp2_names = dir([jd_file_location,'*.jp2']);


size_arr = [numel(required_lat_arr), numel(required_long_arr)];

output_lat_arr = [];
output_long_arr = [];
output_T_arr = [];
output_start_dtm_arr = NaT(0);
output_end_dtm_arr = NaT(0);
output_ltim_arr = [];
output_jd_arr = [];

%% Process files
fprintf('Finding Diviner temperatures for "%s"...\n', ls_name)
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
    
    jd_filename = 'error';
    for jd_file_idx = 1:length(jd_jp2_names)
        if numel(strfind(jd_jp2_names(jd_file_idx).name, tbol_filename(19:27))) > 0
            jd_filename = jd_jp2_names(jd_file_idx).name;
            jd_filename = jd_filename(1:end-4);
        end
    end
    if strcmp('error', jd_filename)
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
    
    if file_idx > 1
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
        
        lat_bool_arr = ismember(lat_arr, required_lat_arr);
        long_bool_arr = ismember(long_arr, required_long_arr);
        bool_matrix = reshape(lat_bool_arr, [], 1) & reshape(long_bool_arr, 1, []);
        lat_arr = lat_arr(lat_bool_arr);
        long_arr = long_arr(long_bool_arr);
        ltim_matrix = ltim_matrix(bool_matrix);
        ltim_matrix = reshape(ltim_matrix, size_arr);
    end
    ltim_matrix = mod(ltim_matrix, 24);
    
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
    tbol_matrix = imread([tbol_file_location, tbol_filename,'.jp2']);
    tbol_matrix = double(tbol_matrix(bool_matrix));
    tbol_matrix = reshape(tbol_matrix, size_arr);

    tbol_matrix = (tbol_matrix*SCALING_FACTOR)+OFFSET;
    tbol_matrix(tbol_matrix == ((MISSING_CONSTANT*SCALING_FACTOR)+OFFSET)) = NaN;
    tbol_matrix(tbol_matrix <= 0) = NaN; % remove obviously unphysical data points
    tbol_matrix(isnan(ltim_matrix)) = NaN;
    
    
    %% Read jd data
    fileID = fopen([jd_file_location, jd_filename,'.lbl'],'r');
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
    jd_matrix = imread([jd_file_location, jd_filename,'.jp2']);
    jd_matrix = double(jd_matrix(bool_matrix));
    jd_matrix = reshape(jd_matrix, size_arr);
    
    jd_matrix = (jd_matrix*SCALING_FACTOR)+OFFSET;
    jd_matrix(jd_matrix == ((MISSING_CONSTANT*SCALING_FACTOR)+OFFSET)) = NaN;
    jd_matrix(jd_matrix < juliandate(datetime(1900,1,1))) = NaN; % remove obviously unphysical data points
    jd_matrix(isnan(ltim_matrix)) = NaN;
    
    
    %% Process data
    output_start_dtm_arr(end+1:end+sum(~isnan(tbol_matrix(:)))) = START_TIME;
    output_end_dtm_arr(end+1:end+sum(~isnan(tbol_matrix(:)))) = STOP_TIME;
    for px_idx = 1:numel(tbol_matrix)
        if ~isnan(tbol_matrix(px_idx))
            output_lat_arr(end+1) = lat_matrix(px_idx);
            output_long_arr(end+1) = long_matrix(px_idx);
            output_T_arr(end+1) = tbol_matrix(px_idx);
%             output_start_dtm_arr(end+1) = START_TIME;
%             output_end_dtm_arr(end+1) = STOP_TIME;
            output_ltim_arr(end+1) = ltim_matrix(px_idx);
            output_jd_arr(end+1) = jd_matrix(px_idx);
        end
    end
    
    if file_idx > progress_idx*(length(tbol_jp2_names)/progress_limit)
        fprintf('#')
        progress_idx = progress_idx + 1;
%         plot_map(lat_arr, long_arr, Tmax_ltim_matrix);
    end
end

fprintf('\n')
output = struct;
output.crater_data = ls_data;
output.lat_arr = output_lat_arr;
if compressed_environment
    output_compressed_long_arr = NaN(size(output_long_arr));
    for long_idx = 1:numel(output_compressed_long_arr)
        [~,value_idx] = min(abs(output_long_arr(long_idx) - compressed_long_arr));
        output_compressed_long_arr(long_idx) = compressed_long_arr(value_idx);
    end
    output.long_arr = output_compressed_long_arr;
    output.uncompressed_long_arr = output_long_arr;
else
    output.long_arr = output_long_arr;
end
output.T_arr = output_T_arr;
output.start_dtm_arr = output_start_dtm_arr;
output.end_dtm_arr = output_end_dtm_arr;
output.ltim_arr = output_ltim_arr;
output.jd_arr = output_jd_arr;
output.dtm_arr = datetime(output_jd_arr, 'convertfrom', 'juliandate');
output.description = 'All Diviner measurements for landing site area';
output.created = datetime;

%% Save outputs
target_path = create_static_path(sprintf('crater_environments/diviner_temperatures/%s_diviner_temperatures', ls_name));
save(target_path, 'output') % save local copy of output as backup

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