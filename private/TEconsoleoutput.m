function TEconsoleoutput(verbosity, message, loglevel, varargin)

% TECONSOLEOUTPUT prints information on program execution to the command 
% line depending on the level of detail requested by the user (verbosity).
% Each output consists of a message and loglevel (level of detail the
% message is considered to have). If the level of verbosity requested by
% the user corresponds to the loglevel, the message is printed.
%
% Input may be a string or cell array of strings. If message is a cell 
% array, it is assumed that the array contains a table and is printed
% accordingly (an additional heading may be passed as varargin), e.g.:
%
%   >> a = {'x='  '4'; 'y=' '8'; 'z=' '3'};
%   >> TEconsoleoutput('info_minor', a, 2, 'my table')
%    base - line NaN: my table
%        x=     4
%        y=     8
%        z=     3
% 
% Messages are indented according to the depth of the call stack and are 
% printed with info on calling function and line number.
%
%
% http://stackoverflow.com/questions/312378/debug-levels-when-writing-an-application
%
% * INPUT
% verbosity - verbosity (string) of command line outputs requested by the 
%             user (see below)
% message   - message to be printed - may be a string or a cell array of
%             strings; if a cell array is provided it is assumed to be a
%             table and printed as such, in this case an optinal heading 
%             (string) can be passed as a varargin
% loglevel  - level of detail (int 0-4) of the message, defined by the
%             program logic (see below)
% varargin  - optional arguments:
%             heading (string) if message is a cell array of strings, i.e.,
%             a table
%
% Verbosity and corresponding logging levels:
%   'none'         - 0 - no output at all
%   'info_major'   - 1 - major program execution steps (e.g. data 
%                        preprocessing, interaction delay reconstruction, 
%                        TE estimation)
%   'info_minor'   - 2 - minor program execution steps (e.g. call to 
%                        subroutines like TEchannelselect)
%   'debug_coarse' - 3 - coarse debug information
%   'debug_fine'   - 4 - very detailed debug information
%
% 
% Version 1.0 by Patricia Wollstadt
% Frankfurt 2015
%

% get current stack for line info in output
stack = dbstack;
if length(stack) > 1
    stack = stack(2:end);
else % if function is called directly
    stack.file = 'base';
    stack.name = 'base';
    stack.line = nan;
end

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
        fprintf('\n')
        error('TRENTOOL ERROR: Unknown verbosity level!')
end

if loglevel > verbosity 
    return
end


if iscell(message)
    message = cell2txt(message, length(stack));
    
    if ~isempty(varargin)
        message = [varargin{1} message];
    end
end


if loglevel == 1
    fprintf(['\n\n%s - line %d: \n' message '\n'], stack(1).file, stack(1).line);         % concatenation instead of %s is needed to make MATLAB do a new line 
        
elseif loglevel == 2 || loglevel == 3
    fprintf('\n')
    for i=1:length(stack)-1;
        fprintf('   ')
    end
    fprintf(['%s - line %d: ' message], stack(1).file, stack(1).line);
    
elseif loglevel == 4
    for i=1:length(stack)
        fprintf('%s - line %d\n', stack(i).file, stack(i).line);
        for j=1:i
            fprintf('\t')
        end
    end
    fprintf(['msg: ' message '\n']);
    
end


function txt = cell2txt(msg, indent)

[rows, cols] = size(msg);
txt = ['\n'];

for r = 1:rows
    
    % add indent for wach row
    for i=1:indent
        txt = [txt '   '];
    end
    
    % add content in each column
    for c = 1:cols
        txt = [txt '\t' msg{r,c}];
    end
    
    % new line for next column
    txt = [txt '\n'];
end

% remove last newline
txt = txt(1:end-2);

