function ret = TEsrmc(pause_time, cmd, varargin)

switch cmd

    case 'maxmem'

        command = 'srmc max gpumem';
        [status,cmdout] = system(command);
        if status == 0    
            gpu_memsize = str2double(cmdout(strfind(cmdout, ' '):end));
            fprintf('Max. GPU memory is %d MB\n', gpu_memsize);
        else
            error('TRENTOOL ERROR: call to srmc returned non-zero exit!')
        end
        
        ret = gpu_memsize;

    case 'request'
        
        resources = varargin{1};

        fprintf(1, '---------------------------------------------------\n');
        fprintf(1, '---------------------------------------------------\n');
        fprintf(1, '%s  requesting %s\n', datestr(now), resources);

        command=sprintf('srmc request %s', resources);

        while (true)
          [status,cmdout] = system(command);
          if (status==0) unit=strtrim(cmdout); break; end
          pause(pause_time+randi(pause_time));
        end

        fprintf(1, '%s  SUCCESS: got unit ''%s''\n', datestr(now), unit);
        fprintf(1, '===================================================\n');
        fprintf(1, '===================================================\n');
        
        if ~strcmp(unit(1:11), '/dev/nvidia')
            error('TRENTOOL error: resource manager returned unkown unit.');
        else
            gpuid = str2double(unit(12:end));
        end
        
        ret = gpuid;

    case 'return'
        
        resources = varargin{1};
        unit      = varargin{2};

        fprintf(1, '---------------------------------------------------\n');
        fprintf(1, '---------------------------------------------------\n');
        fprintf(1, '%s  returning %s\n', datestr(now), resources);

        command=sprintf('srmc return %s %s', unit, resources);

        while (true)
          [status,cmdout] = system(command);
          if (status==0) break; end
          pause(pause_time+randi(pause_time));
        end

        fprintf(1, '%s  resources returned\n', datestr(now));
        fprintf(1, '===================================================\n');
        fprintf(1, '===================================================\n');
    
    otherwise
        
        error('TRENTOOL ERROR: Unkown command for srmc.')
end

