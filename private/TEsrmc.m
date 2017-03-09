function ret = TEsrmc(pause_time, cmd, verbosity, varargin)

%% define logging levels

LOG_INFO_MAJOR = 1;
LOG_INFO_MINOR = 2;
LOG_DEBUG_COARSE = 3;

%%
switch cmd

    case 'maxmem'

        command = 'srmc max gpumem';
        [status,cmdout] = system(command);
        if status == 0
            gpu_memsize = str2double(cmdout(strfind(cmdout, ' '):end));
            msg = sprintf('Max. GPU memory is %d MB', gpu_memsize);
            TEconsoleoutput(verbosity, msg, 2);
        else
            fprintf('\n')
            error('TRENTOOL ERROR: call to srmc returned non-zero exit!')
        end

        ret = gpu_memsize;

    case 'request'

        resources = varargin{1};

        %fprintf(1, '---------------------------------------------------\n');
        %fprintf(1, '---------------------------------------------------\n');
        if ~strcmp(verbosity, 'none')
            TEwaitbar('init',50); fprintf('\b')
            TEwaitbar('init',50); fprintf('\b')
        end
        msg = sprintf('%s  requesting %s', datestr(now), resources);
        TEconsoleoutput(verbosity, msg, 2);

        command=sprintf('srmc request %s', resources);

        while (true)
          [status,cmdout] = system(command);
          if (status==0); unit=strtrim(cmdout); break; end
          pause(pause_time+randi(pause_time));
        end

        msg = sprintf('%s  SUCCESS: got unit ''%s''', datestr(now), unit);
        TEconsoleoutput(verbosity, msg, 2);
        %fprintf(1, '===================================================\n');
        %fprintf(1, '===================================================\n');
        if ~strcmp(verbosity, 'none')
            TEwaitbar('init',50); fprintf('\b')
            TEwaitbar('init',50); fprintf('\b')
        end

        if ~strcmp(unit(1:11), '/dev/nvidia')
            fprintf('\n')
            error('TRENTOOL error: resource manager returned unkown unit.');
        else
            gpuid = str2double(unit(12:end));
        end

        ret = gpuid;

    case 'return'

        resources = varargin{1};
        unit      = varargin{2};

        %fprintf(1, '---------------------------------------------------\n');
        %fprintf(1, '---------------------------------------------------\n');
        if ~strcmp(verbosity, 'none')
            TEwaitbar('init',50); fprintf('\b')
            TEwaitbar('init',50); fprintf('\b')
        end
        msg = sprintf('%s  returning %s', datestr(now), resources);
        TEconsoleoutput(verbosity, msg, 2);

        command = sprintf('srmc return %s %s', unit, resources);

        while (true)
          [status, cmdout] = system(command);
          if (status==0); break; end
          pause(pause_time+randi(pause_time));
        end

        msg = sprintf('%s  resources returned', datestr(now));
        TEconsoleoutput(verbosity, msg, 2);
        %fprintf(1, '===================================================\n');
        %fprintf(1, '===================================================\n');
        if ~strcmp(verbosity, 'none')
            TEwaitbar('init',50); fprintf('\b')
            TEwaitbar('init',50); fprintf('\b')
        end

    otherwise
        fprintf('\n')
        error('TRENTOOL ERROR: Unkown command for srmc.')
end

