function copy_files_to_remote(target_name)
% COPY_FILES_TO_REMOTE copies *.m files from git/code for use on remote
% computers
%
% COPY_FILES_TO_REMOTE copies *.m files from git/code to all remote
% computer folders (remote/pc1, remote/pc2 etc.) if called from git/code,
% or copies *.m files from git/code to current computer's folder otherwise
% (e.g. if called on pc3, copies files to remote/pc3)
%
% COPY_FILES_TO_REMOTE(target_name) copies *.m files from git/code to
% computer with target name (e.g. 'pc1')

if nargin == 0
    computer_name = read_file('computer_name/computer_name.txt');
    computer_name = computer_name{1};
    if strcmpi(computer_name, 'oliverspc2')
        copy_files_to_remote('laptop')
        copy_files_to_remote('pc1')
        copy_files_to_remote('pc2')
        copy_files_to_remote('pc3')
        copy_files_to_remote('pc4')
        return
    else
        target_name = computer_name;
    end
end

if isnumeric(target_name)
    if target_name == 0
        target_name = 'laptop';
    else
        target_name = sprintf('pc%i', target_name);
    end
end

check_dropbox_running

fprintf('Copying files to %s... \t', target_name)
output_folder = sprintf('../../remote/%s', target_name);

original_folder = pwd;
source_folder = '../../git/code';
cd(source_folder)

input_files = dir(fullfile('*.m'));
file_names = { input_files.name };
for k = 1 : length(input_files )
	this_file_name = file_names{k};
	input_file_name = fullfile(pwd, this_file_name);
	output_file_name = fullfile(output_folder, this_file_name);
	copyfile(input_file_name, output_file_name);
end

cd(original_folder)
fprintf('DONE\n')


    function file_text = read_file(code_path)
        % READ_FILE reads in a text file into a cell array with line of the
        % text file as a line in the text array
        fid = fopen(code_path, 'r');
        file_text = {};
        while ~feof(fid)
            file_text{end+1} = fgetl(fid);
        end
        fclose(fid);
    end


    function check_dropbox_running
        % CHECK_DROPBOX_RUNNING checks if dropbox.exe is running, and
        % starts dropbox if it is not already running
        [~,result] = system('tasklist /FI "imagename eq dropbox.exe" /fo table /nh');
        if contains(lower(result), lower('No tasks are running which match the specified criteria'))
            error('Dropbox not running')
        end
    end
end