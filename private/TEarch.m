function [path2TSTOOL] = TEarch

% FUNCTION TEARCH
%
% Adds the folder with TSTool functions to the MATLAB search path. The
% functions checks which OS is running and whether architecture is 32 or 64
% bit.
%
% PW 10/11/2014

% CHANGELOG
% 2015/02/12 - PW: added support for newer MAC versions (.mexmaci64)

% check if TSTool functions are already on the path
h = which('nn_prepare.m');
if ~isempty(h)
    return
end

path2TRENTOOL = which('TEprepare');
path2TRENTOOL = path2TRENTOOL(1:end-11);
path2TSTOOL = fullfile(path2TRENTOOL,'tstool_functions');

extension = mexext;		% gets the mex extension expected by the current architecture

switch extension
  case 'mexw64'
    path2TSTOOL_mex = fullfile(path2TSTOOL,'mex_win', 'mex64');
  case 'mexw32'
    path2TSTOOL_mex = fullfile(path2TSTOOL,'mex_win', 'mex32');
  case 'mexmaci64'
    path2TSTOOL_mex = fullfile(path2TSTOOL,'mex_mac', 'mex64');
  case 'mexmaci32'
    path2TSTOOL_mex = fullfile(path2TSTOOL,'mex_mac', 'mex32');
    case 'mexa64'
    path2TSTOOL_mex = fullfile(path2TSTOOL,'mex_linux', 'mex64');
  case 'mexglx' 
    path2TSTOOL_mex = fullfile(path2TSTOOL,'mex_linux', 'mex32');
  otherwise
    error(['TRENTOOL ERROR: System not supported! Expected extension: ''' mexext '''']);
end

addpath(path2TSTOOL);       % this folder contains the help files for TSTool mex functions
addpath(path2TSTOOL_mex);   % this folder contains the mex functions themselves
