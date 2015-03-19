function [path2TSTOOL] = TEarch

% FUNCTION TEARCH
%
% Adds the folder with TSTool functions to the MATLAB search path. The
% functions checks which OS is running and whether architecture is 32 or 64
% bit.
%
% PW 10/11/2014

% check if TSTool functions are already on the path
h = which('nn_prepare.m');
if ~isempty(h)
    return
end

arch = computer('arch');
path2TRENTOOL = which('TEprepare');
path2TRENTOOL = path2TRENTOOL(1:end-11);

% check OS
if strcmp(arch(1:3),'win')    
    path2TSTOOL = fullfile(path2TRENTOOL,'tstool_functions','mex_win');
else 
    path2TSTOOL = fullfile(path2TRENTOOL,'tstool_functions','mex_linux');
end

addpath(path2TSTOOL);   % this folder contains the help files for TSTool mex functions

% check 32 vs. 64 bit
if strcmp(arch(end-1:end), '64')
    path2TSTOOL_mex = fullfile(path2TSTOOL, 'mex64');
else
    path2TSTOOL_mex = fullfile(path2TSTOOL, 'mex32');
end
    

addpath(path2TSTOOL_mex);   % this folder contains the mex functions themselves