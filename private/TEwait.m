function TEwait(LoopNr)
%
% This function is a very fast and memory efficient job status feedback
% to visually indicate running processes (e.g. 'for' loops)
% It displays a turningbar in brackets at the command prompt. The turning 
% bar works only properly in loops with stepwidth=1.
%
% Input: 
%       LoopNr    = actual loop number
%
% Michael Lindner & Frederic Roux
% Frankfurt am Main, Germany
% 2009

% latest code revision 20/11/09

if LoopNr >= 2; fprintf('\b\b\b'); end
if     mod(LoopNr,4) == 0; rl = '-';
elseif mod(LoopNr,4) == 1; rl = '\';
elseif mod(LoopNr,4) == 2; rl = '|';
else   rl = '/';
end
fprintf('[%s]',rl)
return