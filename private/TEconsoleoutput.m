function TEconsoleoutput(verbosity, message, stack, loglevel)

% Suggestions for logging levels:
% http://stackoverflow.com/questions/312378/debug-levels-when-writing-an-application
%
% log levels:
%   1 - major program execution steps (e.g. data preprocessing, interaction
%       delay reconstruction, TE estimation)
%   2 - minor program execution steps (e.g. call to subroutines like
%       TEchannelselect)
%   3 - coarse debug information
%   4 - fine debug information

switch verbosity
    case 'none'
        verbosity = 0;
    case 'info_major'
        verbosity = 1;
    case 'info_minor'
        verbosity = 2;
    case 'debug_coarse'
        verbosity = 3;
    case 'debug_fine'
        verbosity = 4;
    otherwise
        error('TRENTOOL ERROR: Unknown verbosity level!')
end

if loglevel > verbosity 
    return
end

if loglevel == 1
    fprintf('\n\n%s\n', message);
        
elseif loglevel == 2 || loglevel == 3
    fprintf('\n%s - line %d: %s', stack(1).file, stack(1).line, message);
    
elseif loglevel == 4
    for i=1:length(stack)
        fprintf('%s - line %d\n', stack(i).file, stack(i).line);
        for j=1:i
            fprintf('\t')
        end
    end
    fprintf('msg: %s\n', message);
    
end