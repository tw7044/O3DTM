function crater_data = load_crater_environment(environment_name)
% LOAD_CRATER_ENVIRONMENT loads a crater environment file
%
% crater_data = LOAD_CRATER_ENVIRONMENT(environment_name) loads crater
% environment with name environment_name and returns the data file for the
% environment

crater_data = load(create_static_path(sprintf('crater_environments/%s.mat', environment_name)));
crater_data = crater_data.data;
end