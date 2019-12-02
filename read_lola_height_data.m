function [elevation_matrix, lat_arr, long_arr, mean_elevation] = read_lola_height_data(ppd, lat_limits, long_limits)
% READ_LOLA_HEIGHT_DATA loads LOLA elevation data
%
% [elevation_matrix, lat_arr, long_arr, mean_elevation] = READ_LOLA_HEIGHT_DATA(ppd) reads and
% returns LOLA height for given ppd

if nargin == 0 || ppd == 4
    image_location = create_static_path('4_ldem/ldem_4.jp2');
    filelocation = create_static_path('4_ldem/');
    filename = 'ldem_4_jp2';
elseif ppd == 16
    image_location = create_static_path('16_ldem/ldem_16.jp2');
    filelocation = create_static_path('16_ldem/');
    filename = 'ldem_16_jp2';
elseif ppd == 128
    image_location = create_static_path('128_ldem/ldem_128.jp2');
    filelocation = create_static_path('128_ldem/');
    filename = 'ldem_128_jp2';
else
    error('Unexpected ppd')
end
if nargin < 2 || isempty(lat_limits)
    lat_min = -inf;
    lat_max = inf;
else
    lat_min = min(lat_limits);
    lat_max = max(lat_limits);
end
if nargin < 3 || isempty(long_limits)
    long_min = -inf;
    long_max = inf;
else
    long_min = min(long_limits);
    long_max = max(long_limits);
end

fileID = fopen([filelocation,filename,'.lbl'],'r');
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
OFFSET = find_str(text_data,pointer(1));

pointer = strfind(text_data,'SCALING_FACTOR');
SCALING_FACTOR = find_str(text_data,pointer(1));


% fprintf('\n\nfile labels: \n MAP_RESOLUTION \t: %d\n CENTER_LONGITUDE \t: %d\n CENTER_LATITUDE \t: %d\n',MAP_RESOLUTION,CENTER_LONGITUDE,CENTER_LATITUDE)
% fprintf(' LINE_PROJECTION_OFFSET \t: %d\n SAMPLE_PROJECTION_OFFSET \t: %d\n\n MAXIMUM_LATITUDE \t: %d\n',LINE_PROJECTION_OFFSET,SAMPLE_PROJECTION_OFFSET,MAXIMUM_LATITUDE)
% fprintf(' MINIMUM_LATITUDE \t: %d\n WESTERNMOST_LONGITUDE \t: %d\n EASTERNMOST_LONGITUDE \t: %d\n OFFSET \t: %d\n SCALING_FACTOR \t: %d\n',MINIMUM_LATITUDE,WESTERNMOST_LONGITUDE,EASTERNMOST_LONGITUDE,OFFSET,SCALING_FACTOR)
fclose(fileID);


elevation_matrix = imread(image_location);
lat_arr = NaN(1,size(elevation_matrix,1));
long_arr = NaN(1,size(elevation_matrix,2));
for LINE = 1:size(elevation_matrix,1)
    lat_arr(LINE) = CENTER_LATITUDE - (LINE - LINE_PROJECTION_OFFSET - 1) / MAP_RESOLUTION;
end
for SAMPLE = 1:size(elevation_matrix,2)
    long_arr(SAMPLE) = CENTER_LONGITUDE + (SAMPLE -SAMPLE_PROJECTION_OFFSET - 1) / MAP_RESOLUTION;
end

%% Get region data
lat_bool_arr = (lat_arr >= lat_min) & (lat_arr <= lat_max);
long_bool_arr = (long_arr >= long_min) & (long_arr <= long_max);
elevation_bool_matrix = reshape(lat_bool_arr, [], 1) & reshape(long_bool_arr, 1, []);
lat_arr = lat_arr(lat_bool_arr);
long_arr = long_arr(long_bool_arr);
elevation_matrix = elevation_matrix(elevation_bool_matrix);
elevation_matrix = reshape(elevation_matrix, numel(lat_arr), numel(long_arr));
elevation_matrix = double(elevation_matrix);
elevation_matrix = (elevation_matrix*SCALING_FACTOR) + OFFSET;
mean_elevation = mean(elevation_matrix(:));
end

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