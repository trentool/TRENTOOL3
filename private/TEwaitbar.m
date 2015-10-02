function TEwaitbar(varargin)

% TEWAITBAR: Plots a waitbar ('-----') with the same indentation as 
% TEconsoleoutput.
%
% * INPUT PARAMETERS 
%
%   TEwaitbar('init', wait_length, verbosity)
%        initialise a reference waitbar with length wait_length (bar is not 
%        plotted if verbosity is set to 'none')
%
%   TEwaitbar('update', iteration, verbosity)
%        update wait bar by plotting one dash ('-')
%
% Version 1.0 by Patricia Wollstadt, Frankfurt 2015


status = varargin{1};
if strcmp(varargin{end}, 'none')
    return
end

stack = dbstack;

switch status
    
    case 'init'
        
        wait_length = varargin{2};
        
        fprintf('\n')
        for i=1:length(stack)-2;
            fprintf('   ')
        end
        for ii=1:wait_length
            fprintf('-')
        end
        fprintf('\n')
        
    case 'update'
        
        iteration = varargin{2};
        
        if iteration == 1
        for i=1:length(stack)-2;
            fprintf('   ')
        end
        end
        fprintf('-');
end