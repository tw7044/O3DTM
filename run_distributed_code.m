function run_distributed_code(code_path, always_on)
% RUN_DISTRIBUTED_CODE runs code specified in code_path on multiple
% machines, using Dropbox to ensure different machines coordinate which
% lines of code they execute.
% Specify distributed code in a .txt file at code_path, each line of the
% file can contain a series of MATLAB statements that are executed on a
% single machine (in series). Each seperate line of the .txt file may be
% executed on seperate machines (in parallel). The seperate lines of code
% should therefore be independent of each other and be able to execute out
% of sequence. Lines of code in the .txt file should be given in order of
% precedence as earlier lines will be executed first.
% During execution, this script will edit the .txt file, adding tags (e.g.
% $RUNNING{...}) to each line of code to indicate if it should be run by
% other machines. A specific computer can be prevented from running a line
% of code by adding $AVOID{computer_name} to the line of code that computer
% should not run.
% After execution, this script will add either $DONE{computer_name} or
% $ERROR{computer_name} to the line of code to indicate if the line was
% executed sucessfully. Once all lines are marked as $DONE{...} (or a line
% has encountered an error on all machines), the distributed execution is
% complete.
% Additional lines of code can be added to the END of the .txt file during
% execution of this script if necessary.
%
% RUN_DISTRIBUTED_CODE runs code located in default file
%
% RUN_DISTRIBUTED_CODE(code_path) defines code file location
% 
% RUN_DISTRIBUTED_CODE(code_path, always_on) defines if script should wait
% for new commands once all distributed code is run

if nargin < 1
    code_path = '../../remote/distributed_code.txt';
end
if nargin < 2
    always_on = true;
end
computer_name = read_file('computer_name/computer_name.txt');
computer_name = computer_name{1};

warning('on', 'all') % ensure any warnings will be shown to user

fprintf('RUNNING DISTRIBUTED CODE\n\tcode_path: \t\t%s\n\tcomputer_name: \t%s\n', code_path, computer_name)

%% Run code in three stages
% First stage - run each line once
% Second stage - rerun lines which have encountered errors
% Third stage - run lines which haven't finished execution on spare
% machines (so that one incredibly slow/frozen machine doesn't hold up
% execution everywhere)
% If always_on is set to true, code will keep running, periodically
% checking if any new tasks are added to the text file
while true
    for stage_idx = 1:3
        fprintf('\n\n')
        fprintf('************************************************\n')
        fprintf('***** RUNNING DISTRIBUTED CODE (STAGE %i/3) *****\n', stage_idx)
        fprintf('************************************************\n')
        continue_execution = true;
        while continue_execution
            % Process the code file and execute code contained in it
            continue_execution = process_file(code_path, computer_name, stage_idx);
        end
    end
    fprintf('\n\n')
    fprintf('************************************************\n')
    fprintf('******** DISTRIBUTED EXECUTION COMPLETE ********\n')
    fprintf('************************************************\n\n')
    if always_on
        fprintf('******************************\n')
        fprintf('** Waiting for new commands **\n')
        for time_idx = 1:30
            fprintf('*')
            pause(60*2) % check for new code every hour
        end
        fprintf('\n')
    else
        break
    end
end



%% Define functions
    function continue_execution = process_file(code_path, computer_name, stage_idx)
        % PROCESS_FILE reads in the code file and processes it line-by-line
        % to find a line of code that needs executing
        
        file_text = read_file(code_path, computer_name);
        continue_execution = false;
        
        for line_idx = 1:numel(file_text)
            % Progressively test each line of the code file to see if it is
            % currently executing on another machine, and if not, execute
            % it here
            
            file_line = file_text{line_idx};
            
            if numel(file_line) > 0 && ~strcmp(file_line(1), '%')
                % Only process lines which could contain executeable
                % content
                
                if read_tag(file_line, 'avoid', computer_name)
                    continue % skip over line if set to avoid this computer
                    
                elseif tag_exists(file_line, 'done')
                    continue % skip over line if already done
                    
                elseif stage_idx == 1 && tag_exists(file_line, 'error')
                    continue % skip over line (for now) if has caused error
                    
                elseif stage_idx > 1 && read_tag(file_line, 'error', computer_name)
                    continue % skip over line which has caused erorr on this machine
                    % line will be attempted at stage 2 if it has caused error
                    % elsewhere (as may crash on low memory machine)
                    
                elseif stage_idx < 3 && tag_exists(file_line, 'running') && ~read_tag(file_line, 'error', read_tag(file_line, 'running'))
                    continue % skip over line (for now) if running elsewhere and not
                    % an error line will be attempted at stage 3 if it is
                    % still running elsewhere (as may be stuck on slow
                    % machine)
                    
                else
                    % execute the code contained in this line here
                    
                    % mark the line as executing on this maching
                    file_line = add_tag(file_line, 'running', computer_name);
                    file_text{line_idx} = file_line;
                    write_file(code_path, file_text);
                    
                    % execute the code contained in the file line
                    code_success = execute_code(file_line);
                    
                    % deal with sucessful code execution and errors in the
                    % code execution. Re-read in file text to ensure any
                    % changes during code execution are delt with.
                    file_text = read_file(code_path, computer_name);
                    file_line = file_text{line_idx};
                    
                    if code_success
                        file_line = add_tag(file_line, 'done', computer_name);
                        file_text{line_idx} = file_line;
                        write_file(code_path, file_text);
                    else
                        file_line = add_tag(file_line, 'error', computer_name);
                        file_text{line_idx} = file_line;
                        write_file(code_path, file_text);
                    end
                    
                    % execution will have taken time, so break here so that
                    % the code file is re-read so that any outstanding
                    % lines can be run
                    continue_execution = true;
                    break
                    
                end
            end
        end
    end

    function file_text = read_file(code_path, computer_name)
        % READ_FILE reads in a text file into a cell array with line of the
        % text file as a line in the text array
        check_dropbox_running
        if nargin > 1
            force_consistent_timing(computer_name)
        end
        fid = fopen(code_path, 'r');
        file_text = {};
        while ~feof(fid)
            file_text{end+1} = fgetl(fid);
        end
        fclose(fid);
    end

    function write_file(code_path, file_text)
        % WRITE_FILE writes to a text file from a cell array with each line
        % of the text file from a cell in the array
        fid = fopen(code_path, 'w');
        for line_idx = 1:numel(file_text)
            fprintf(fid, strrep(file_text{line_idx},'%','%%'));
            fprintf(fid,'\r\n');
        end
        fclose(fid);
    end

    function force_consistent_timing(computer_name)
        % FORCE_CONSISTENT_TIMING pauses execution of code until a specific
        % minute (e.g. pc3 has minutes ending in 5). This is used so that
        % each different machine has a different window of time it is
        % allowed to interact with the text file to help avoid clashes
        % where two machines are trying to interact with the file at once.
        % In this sub function, 'time' refers to the final digit of the
        % current minute, e.g. 12:32pm has a 'time' of 2.
        switch computer_name
            case 'AOPPWT4'
                io_time = 1;
                % i.e. pc1 can only interact with the text during a minute
                % long window file at 1, 11, 21, 31, 41 & 51 minutes past
                % every hour.
            case 'AOPPWT3'
                io_time = 2;
            case 'AOPPWT11'
                io_time = 3;
            case 'AOPPWT14'
                io_time = 4;               
            case 'AOPPWT15'
                io_time = 5;               
            case 'AOPPWT16'
                io_time = 6;                               
            case 'AOPPWT10'
                io_time = 7;                                               
            case 'AOPPWT17'
                io_time = 8;                                                               
            case 'mainpc'
                return
            otherwise
                io_time = 9;
        end
        current_dtm = datetime; % get current datetime
        current_time = mod(minute(current_dtm),10); % get current 'time'
        current_sec = second(current_dtm); % get current seconds
        if current_time ~= io_time
            % pause execution until start of next 'time' window if not
            % currently in a window
            pause(mod(60*(io_time-current_time)-current_sec, 10*60));
        end
    end

    function tag_value = read_tag(file_line, tag_name, tag_value_test)
        % READ_TAG tests for the existance of a tag of the form
        % $TAG_NAME{tag_value} in a line in the text file. Returns false if
        % the tag doesn't exist and returns tag_value if the tag does exist
        tag_str = sprintf('$%s{', tag_name);
        tag_idx_arr = strfind(upper(file_line), upper(tag_str));
        if numel(tag_idx_arr) == 0
            tag_value = false;
            return
        end
        if nargin == 3
            for tag_idx = tag_idx_arr
                tag_start = strfind(file_line(tag_idx:end), '{');
                tag_end = strfind(file_line(tag_idx:end), '}');
                tag_value = file_line(tag_idx+tag_start(1):tag_idx+tag_end(1)-2);
                tag_value = strcmpi(tag_value, tag_value_test);
                if tag_value
                    break
                end
            end
        else
            tag_idx = tag_idx_arr(1);
            tag_start = strfind(file_line(tag_idx:end), '{');
            tag_end = strfind(file_line(tag_idx:end), '}');
            tag_value = file_line(tag_idx+tag_start(1):tag_idx+tag_end(1)-2);
        end
    end

    function result = tag_exists(file_line, tag_name)
        % TAG_EXISTS returns boolean to indicate if a tag of the form
        % $TAG_NAME{tag_value} exists in the file_line
        if read_tag(file_line, tag_name)
            result = true;
        else
            result = false;
        end
    end

    function file_line = add_tag(file_line, tag_name, tag_value)
        % ADD_TAG adds tag of the form $TAG_NAME{tag_value} to the end of a
        % line of text
        file_line = sprintf('%s $%s{%s}', file_line, upper(tag_name), lower(tag_value));
    end

    function success = execute_code(file_line)
        % EXECUTE_CODE extracts matlab code from the start of a a line in
        % the text file and then executes the code. Returns boolean
        % indicating if the code execution was sucessful (true) or
        % contained any errors (false)
        code_text = strsplit(file_line, '$');
        code_text = code_text{1};
        success = true;
        rehash
        fprintf('\n> %s\n', code_text)
        try
            eval(code_text)
        catch me
            warning(getReport(me, 'extended', 'hyperlinks', 'on'))
            success = false;
        end
    end

    function check_dropbox_running
        % CHECK_DROPBOX_RUNNING checks if dropbox.exe is running, and
        % starts dropbox if it is not already running
        [~,result] = system('tasklist /FI "imagename eq dropbox.exe" /fo table /nh');
        if contains(lower(result), lower('No tasks are running which match the specified criteria'))
            winopen('C:\Program Files (x86)\Dropbox\Client\Dropbox.exe') % launch dropbox
            fprintf('\n*** Starting Dropbox ... ***\n')
            pause(60*5) % give dropbox time to get running and sync
        end
    end
end