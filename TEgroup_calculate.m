function TEgroup_calculate(filename)

% TEGROUP_CALCULATE: This function calculates the transfer entropy values
% for the 
%
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% !!!                   The function TEgroup_prepare                    !!!
% !!!               has to be run on all datasets first!                !!!
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%
% You can call this function directly as follows:
%         TEgroup_calculate(filename)
%
%
% * REFERENCE INFORMATION
%
%   - transfer entropy
%     - The concept of TE appears in Schreiber's article,
%       "Measuring Information Transfer", Phys. Rev. Lett. 85, 461 - 464
%       (2000).
%     - For the estimation of probability densities needed for the TE
%       computation, the function implements the Kraskov-Stoegbauer-
%       Grassberger estimator described in Kraskov et al. "Estimating
%       mutual information", Phys. Rev. E 69 (6) 066138, (2004).
%
%   - permutation test
%       - Maris & Oostenveld (2007). Nonparametric statistical testing of
%         EEG- and MEG-data. J. of Neuroscience Methods, 164, 177-190.
%
%
% * DEPENDENCIES
%     - Package TSTOOL is used at nearest neighbors searches
%       required for the KSG estimator. (Gnu Public License)
%       http://www.dpi.physik.uni-goettingen.de/tstool/
%     - The following Matlab toolboxes:
%         - statistic toolbox
%     - The functions
%         - transferentropy
%         - TEcmc
%         - TEperm
%         - TEvalues
%         - TEwait
%
% 
%
% OUTPUT PARAMETERS
%
%
%  TEresult             = Output structure of the function tranferentropy
%          .TEmat       = resultmatrix including transfer entropy(TE)
%                         values (channelpairs x u x trial)
%          .MImat       = resultmatrix including mutual information (MI)
%                         values (channelpairs x u x trial)
%          .dimord      = 'channelpair_u_trial'; the dimensions of TEmat 
%                         and MImat
%          .cfg         = configuration file used to calculate TE
%          .trials      = trial numbers selected from raw dataset
%          .act         = ACT matrix (trial x (target) channel)
%          .sgncmb      = labels of channel combinations (source -> target)
%          .TEprepare   = results of the function TEprepare from the
%                         data
%   if instantaneous mixing is found in the data, then another field will 
%   be added:
%          .instantaneousmixing = matrix (channelcombi) which indicates were
%          the instantaneous mixings were found (1) or not (0).%
%
%
%
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation;
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY;
%
% Version 1.0 by Michael Lindner, Raul Vicente, Michael Wibral
% Frankfurt 2010
%

% CHANGELOG
%
% 2013-06-21 PW: corrected data handling in .TEpermvalues for cases of 
% volume conduction
% 2013-06-21 PW: added the hanndling of the field 'cfg.extracond'
% 2013-06-25 PW: individual embedding delays are now used for all subjects


warning('!!!!!!!! THIS FUNCTION IS DEPRECATED - USE TEgroup_prepare and TEgroup_stats instead !!!!!!!!');
pause(5);



%% Remember the working directory
working_directory1 = pwd;


%% load data
% -------------------------------------------------------------------------
fprintf('\nLoad data');

varinfile = who('-file',filename);
load(filename)
x = strcat('data=',varinfile{1},';');
eval(x)    

fprintf(' - ok');



%% check data
% -------------------------------------------------------------------------
fprintf('\nCheck data');

if ~isfield(data, 'TEprepare') 
     error('TRENTOOL ERROR: The functions TEprepare and TEgroup_prepare must be performed on all datasets first, see help!');
end

if ~isfield(data, 'TEgroupprepare') 
     error('TRENTOOL ERROR: The function TEgroup_prepare must be performed on all datasets first, see help!');
end

fprintf(' - ok');


%% get cfg from TEgroupprepare 
% -------------------------------------------------------------------------



cfg = data.TEgroupprepare.cfg;

% use individual tau for TE calculation
if strcmp(data.TEprepare.cfg.optimizemethod, 'ragwitz');
    cfg.tau = data.TEprepare.opttau;
end
    
% handle cfg field 'extracond'
if ~isfield(data.TEgroupprepare.cfg,'extracond')
    cfg.extracond = 'no';
end


%% calculate TE
% -------------------------------------------------------------------------

cfg.calctime = 'yes';

fprintf('\nStart calculating transfer entropy for unshuffled data');
cfg.shuffle = 'no';
[TEresult] = transferentropy(cfg,data);
TEresult.TEprepare = data.TEprepare;




%% perform shift test
% -------------------------------------------------------------------------
% TEshift is created inside transferentropy.m as a reduced version of
% TEresult without certain fields. TEshift is never written to disk/file
% to avoid later confusion. Please save TEshift yourself if necessary.

if strcmp(cfg.shifttest, 'yes')
    
    % calculate TE for shifted data
    fprintf('\nStart calculating transfer entropy for shifted data');    
    cfg.calctime = 'no';
    cfg.numpermutation = 190100;
    cfg.shuffle = 'no';
    [TEshift] = transferentropy(cfg,data,'shifttest');
    
    % set parameters for shift test
    fprintf('\nStart permutation tests for shift test');
    cfg.permstatstype = 'indepsamplesT';
    cfg.tail = 1;
    cfg.correctm = 'FDR';
    
    if strcmp(cfg.shifttesttype, 'TE>TEshift')
        
        % do permutation test
        cfg.alpha   = 0.05;        
        TEpermshift = TEperm(cfg,TEresult,TEshift);
        
        % evaluate results from permutation test
        fprintf('\nanalyze shift test\n');
        indexinstmix = find(TEpermshift.TEpermvalues(:,2)==0);
        if size(indexinstmix,1) == 0
            fprintf('No instantaneous mixing found!\n')
        else
            fprintf(strcat(num2str(size(indexinstmix,1)),' instantaneous mixings found by strict shifttest!\nFor these cases TEvalues of all trials are set to NaN!\n'))
            mask=repmat((TEpermshift.TEpermvalues(:,2)-1)*-1, [1 1 size(TEresult.TEmat,2)]);
            TEresult.TEmat(mask==1) = NaN;
            TEresult.MImat(mask==1) = NaN;
            TEresult.instantaneousmixing = (TEpermshift.TEpermvalues(:,2)-1)*-1;
        end
        
        
    elseif strcmp(cfg.shifttesttype, 'TEshift>TE')
        
        % do permutation test
        cfg.alpha   = 0.05/length(data.TEgroupprepare.files);        
        TEpermshift = TEperm(cfg,TEshift,TEresult);
        
        % evaluate results from permutation test
        fprintf('\nanalyze shift test\n');
        indexinstmix = find(TEpermshift.TEpermvalues(:,2)==1);
        if size(indexinstmix,1) == 0
            fprintf('No instantaneous mixing found!\n')
        else
            fprintf(strcat(num2str(size(indexinstmix,1)),' instantaneous mixings found by non-strict shifttest!\nFor these cases TEvalues of all trials are set to NaN!\n'))
            mask=repmat(TEpermshift.TEpermvalues(:,2), [1 1 size(TEresult.TEmat,2)]);
            TEresult.TEmat(mask==1) = NaN;
            TEresult.MImat(mask==1) = NaN;
            TEresult.instantaneousmixing = TEpermshift.TEpermvalues(:,2);
        end
    end  
    
    clear TEpermshift mask
    
    cfg = rmfield(cfg, 'alpha');
    cfg = rmfield(cfg, 'correctm');
    cfg = rmfield(cfg, 'permstatstype');
    cfg = rmfield(cfg, 'numpermutation');
    cfg = rmfield(cfg, 'tail');
    cfg = rmfield(cfg, 'calctime');
end

%% append .TEgroupprepare, .TEprepare, .cfg
% -------------------------------------------------------------------------

TEresult.TEgroupprepare = data.TEgroupprepare;
TEresult.TEprepare      = data.TEprepare;
TEresult.cfg            = cfg;

%% save data
% -------------------------------------------------------------------------

savename = strcat(filename(1:end-25),'TE_output.mat');
save(savename, 'TEresult');

%% Returning to the working directory
cd(working_directory1)











