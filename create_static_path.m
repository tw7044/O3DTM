function static_path = create_static_path(internal_path)
% CREATE_STATIC_PATH returns location where large data files are stored
%
% CREATE_STATC_PATH(internal_path) returns file location with
% static_path/internal_path
%
% CREATE_STATIC_PATH returns folder locatio of static folder


if nargin == 0
    internal_path = '';
end
static_root = '../../mphys_project_static/';
static_path = [static_root, internal_path];
end