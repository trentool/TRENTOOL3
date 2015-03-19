function [ TEbinom ] = TEsurrogate_binomstats( cfg, FilesCell)
%TESURROGATE_BINOMSTATS: Calculates the binomial statistic for the presence
% of a link over a set of results or results files from a preceeding calls
% to TEsurrogatestats
%
% * DEPENDENCIES
%     - The following Matlab toolboxes:
%         - statistics toolbox
%
% You can call this function directly as follows:
%       [ TEbinom ] = TEsurrogate_binomstats(cfg, FilesCell)  
%
% * INPUT PARAMETERS
%
%   cfg: The configuration CAN contain:
%
%   cfg.alpha   = statistical signifiance threshold (default = 0.05)
%
%   FilesCell: The data can be passed in two ways:
%
%   1. a cell array of structures of dimension Nx1 or 1xN, where N is the number
%      of elements in the set, eachstructure in a cell MUST contain the
%      following fields (e.g. you can put in each cell the output 
%      obtained from a call to InteractionDelayReconstruction_analyze.m or 
%      a call to TEsurrogatestats.m):
%
%       .sgncmb       = sgncmb signal combination 
%
%       .TEpermvalues = matrix with size channelpair x 6
%                           The second dimension includes (row-wise):
%                           1 - p_values of the statistic within the
%                               distribution given by the permutations
%                           2 - 1 (0), if the statistics is significant at
%                               the prescribed alpha level (or not)
%                           3 - 1 (0), if the statistics is significant
%                               after correction for mulitple comparisons
%                               (or not)
%                           4 - 1 (0), mean difference or tvalue of mean
%                               difference depending on cfg.permstatstype
%                           5 - 1 (0), if instantaneous mixing (volume
%                               conduction) exists (or not)
%                           6 - delay times u (not mandatory)
%       .cfg         = configure used to compute the TE, the folowing
%                       fields are needed:                       
%                           .alpha statistical signifiance level
%
%   2. a cell array of dimensions Nx1 or 1xN, where N is the number of elements 
%       in the set, where each cell includes the name of a file of the set
%
% * OUTPUT PARAMETERS 
%   TEbinom
%       .TEpermvalues = matrix with size signalombinations x 6 (for the exact
%                       specification see INPUT PARAMETERS), with the
%                       following fields:
%                           1 - p-value from the binomial test
%                           2 - 1 (0), if the statistics is significant at
%                               the prescribed alpha level (or not)
%                           3 - 1 (0), if the statistics is significant
%                               after correction for mulitple comparisons
%                               (or not)
%                           4 - mean TE value for non NaN data
%                           5 - 1 (0), if instantaneous mixing (volume
%                               conduction) exists (or not) over all of the
%                               elements, on the same channel
%                           6 - medium delay times, over all of the elements
%       .sgncmb    = signal combination
%       .TEsteps   = string field containing the name of the functions used
%                    in the workflow
%       .occurrences   = vector with counts of occurences of individual
%                        links in all subjects
%       .bino_cdf      = return binomical cdf for all values from 1 to 
%                        total number of cases
%       .occurr_thresh = threshold for no. occurrences that make a link
%                        statistically significant given the alpha and the
%                        total number of cases
%
% NP - 15/06/2012
%
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation;
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY;
%
% Version 1.0 by Nicu Pampu,Michael Lindner, Raul Vicente, Michael Wibral
% Frankfurt 2012



% CHANGELOG
%
% 2012-06-15: NP creating file
% 2012-07-04: NP adding documentation
% 2014-04-04: MW, bugfix - significance was corrected by number of subjects
% instead of number of channels
% 2014-10-21: PW added three output parameters
% 2014-10-21: PW, bugfix - alpha from cfg was not used
% 2015-01-19: PW, bugfix - 0 occurrences caused function to crash
% 2015-01-29: PW, occurrences thresholds needed a - 1 (the bino cdf starts with 0 occurences)

if ~isfield(cfg,'alpha')
    sign_level = 0.05;      %the default value for binomial statistic signifiance
else 
    sign_level = cfg.alpha;
end

%% check and prepare data 
data = prepare_data(FilesCell); %prepare the data file for further use

binomial_q = data{1,1}.cfg.alpha; %set signifiance level 

%% process the binomial teston group data

count_occurence = zeros(size(data{1}.TEpermvalues,1),1);
nr_nan   = zeros(size(data{1}.TEpermvalues,1),1);
count_TE =zeros(size(data{1}.TEpermvalues,1),1);

if size(data{1}.TEpermvalues,2) == 6
    count_u  = zeros(size(data{1}.TEpermvalues,1),1);
end

for i = 1:size(data,2)
    for j = 1:size(data{1}.TEpermvalues,1)         
        if data{1,i}.TEpermvalues(j,2)
            count_occurence(j) = count_occurence(j) +1; %occurence of signifiance but not corected over the subjects
        end
        if (data{1,i}.TEpermvalues(j,5)>0) || (data{1,i}.TEpermvalues(j,5)<0)
            nr_nan(j) = nr_nan(j)+1;  %occurence of nan's over the subjects 
        else 
            count_TE(j) = count_TE(j) + data{1,i}.TEpermvalues(j,4);     
        end
        if size(data{1}.TEpermvalues,2) == 6
            count_u(j) = count_u(j)+data{1,i}.TEpermvalues(j,6);
        end
    end
    
end

numcmb = size(data{1}.TEpermvalues,1); % MW 2014-04-04 - for bugfix below
total_occ = size(data,2);
TEbinom.TEpermvalues = zeros(size(data{1}.TEpermvalues,1),6);

for i = 1:size(data{1}.TEpermvalues,1)
        TEbinom.TEpermvalues(i,1) =  binomtest(count_occurence(i),(total_occ-nr_nan(i)),binomial_q);
        if TEbinom.TEpermvalues(i,1) <= sign_level%signifiance uncorected
            TEbinom.TEpermvalues(i,2) = 1;
        end
        if TEbinom.TEpermvalues(i,1) <= sign_level/numcmb; % MW 2014-04-04 - bugfix
            TEbinom.TEpermvalues(i,3) = 1; %signifiance corrected
        end
        TEbinom.TEpermvalues(i,4) = count_TE(i)/(total_occ-nr_nan(i));    %Mean TE over the subjects
        if (total_occ-nr_nan(i)) == 0
            TEbinom.TEpermvalues(i,5) = 1;      %Volume Conduction field
        end
        if size(data{1}.TEpermvalues,2) == 6
            TEbinom.TEpermvalues(i,6) = (count_u(i))/(total_occ-nr_nan(i));             %mean U value - 1 (exclude 0 value)
        end
end

TEbinom.cfg = cfg;                      %store configuration
TEbinom.sgncmb = data{1,1}.sgncmb;      %store sgncmb for further use
TEbinom.TEpermvaluesTmp = data;         %store old data
TEbinom.occurrences =  count_occurence; %return vector with counted occurences            
TEbinom.bino_cdf = 1-binocdf(1:total_occ,total_occ,binomial_q);     % return binomical cdf for all values from 1 to total number of cases
TEbinom.occurr_thresh = find(TEbinom.bino_cdf < binomial_q,1) - 1;  % no. occurrences needed to make a link significant

if ~isfield(data{1,1},'TEsteps')      %adding structure with changings
    TEbinom.TEsteps = 'TEsBs';
else TEbinom.TEsteps = strcat(data{1,1}.TEsteps,'_TEsBs');
end


end
function data = prepare_data(TMPdata)
%PREPARE_DATA is used to chose between the input formats, it passes the
%correct data for the binomial group test 
if ischar(TMPdata{1,1})   %check for data type, if it is char(containing file path) or
    state = 1;
elseif isfield(TMPdata{1,1},'TEpermvalues')           %if it is data   
    state = 2;
else error('TRENTOOL error:wrong data input, see help')
end

if size(TMPdata,1)>1
    TMPdata = TMPdata';
    if size(TMPdata,1)>1
        error('TRENTOOL error:data input must be in 1xN or Nx1 format')
    end
end

%% Remember the working directory
working_directory1 = pwd;

if state == 1
% load Data
    DataCell={};
    for ll = 1:length(TMPdata)
        varinfile = who('-file',TMPdata{ll});
        load(TMPdata{ll});
        x = strcat('DataCell{ll}=',varinfile{1},';');
        eval(x)
        y=strcat( ['clear ',varinfile{1} ]);
        eval(y)
        clear x y varinfile
    end
    clear ll

    nrdata = length(DataCell);
    
    if nrdata ~= length(TMPdata)
        error('TRENTOOL error: unequal number of loaded Data and entries in FileCell')
    end

    check_data(DataCell);
    data = DataCell;   
end
if state == 2
    check_data(TMPdata);   
    data = TMPdata;
end
cd(working_directory1);
end
function check_data(data)
%this function checks for the next types of correctness:
%       -test for field TEpermvalues to be the same size for all of the 
%       entries
%       -test for the sgncmb to be the same and present in all of the
%       entries
%       -test for field cfg.alpha to contain the same value in all of the
%       entries

nr_comb = size(data{1,1}.TEpermvalues,1);
binomial_q = data{1,1}.cfg.alpha;

for i = 1:size(data,2)
    if (~isfield(data{1,i},'sgncmb'))||(size(data{1,i}.sgncmb,1) ~= nr_comb)|| (size(data{1,i}.sgncmb,2)~=2)
        error(strcat('TRENTOOL error: There are errors in datafile. Check sgncmb field:',num2str(i)));
    end
end
for i=2: size(data,2)
    aux = size(data{1,i}.TEpermvalues) ~= size(data{1,i-1}.TEpermvalues);
    if aux(1)||aux(2)
         error(strcat('TRENTOOL error: There are errors in datafile. Check TEpermvalues field:',num2str(i)));
    end    
    if data{1,i}.cfg.alpha ~= binomial_q
        error(strcat('TRENTOOL error: There are errors in datafile. Check cfg.alpha'));
    end
end
for i=2: size(data,2)
    for j = 1:nr_comb
        if ~strcmp(data{1,i}.sgncmb(j,1),data{1,i-1}.sgncmb(j,1))||~strcmp(data{1,i}.sgncmb(j,2),data{1,i-1}.sgncmb(j,2))
         error(strcat('TRENTOOL error: Mismatch sgncmb label in entry:',num2str(i)));   
        end
    end
end

end
function [ bprob ] = binomtest(occur, tot_nr, alpha)
% Performs a binomial test to determine the probablility of outcome given
% the number of succesful outcomes, total number of outcomes and a prop
% Use as
%    [bprob] = binomtest(occur, tot_nr, alpha)
%
% Input:
%       occur - The observed numebr of successful outcomes
%       tot_nr- The total number of outcomes (successful or not)
%       alpha - The hypothesis probability of succes.
% Output:
%       prob  - The probability of getting this and more of this number of
%       outcomes from the totatl number of outcomes giving the probability
%       of succes
%  Dependences: Matlab stats toolbox

    if occur == 0
        bprob=1-binocdf(0,tot_nr,alpha);
    else
        bprob=1-binocdf(occur-1,tot_nr,alpha); % occur-1 for >=, binocdf is accepting only >
    end
end
