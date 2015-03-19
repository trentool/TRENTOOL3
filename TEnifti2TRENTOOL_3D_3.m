function TEnifti2TRENTOOL_3D_3(cfg)

% TENIFTI2TRENTOOL - This functions converts the fMRI data of ROIs from SPM
% nifti format into the rawdata format of Fieldtrip to use it in TRENTOOL.
% This function creates a dataset for each subject using eather the time 
% series of a set of 3x3x3 volumes as single trials or the single time
% series of the voxel and surrounding 8 voxels for embedding. 
% Therefore a mask with one voxel
% or the individual peak of a mask can be used. This function is designed
% to use the outputs of SPM as input for masks and contrast images. This
% function is not tested for nifti files of other data analyses software
% packages.
%
% Use the function as followed:
%           TENifti2TRENTOOL_3D(cfg)
%
%
% REQUIREMENT:
%   fMRI data has to be organized in the following way:
%
%   |_folder of subject data
%   |   |_folder subject 1
%   |   |   |_run1     <- names of subdir must be the same for all subjects
%   |   |   |_run2        in case of only one run: no subdirectory!
%   |   |   ...
%   |   |_folder subject 2
%   |   |   |_run1
%   |   |   |_run2
%   |   |   ...
%   |   ...
%   |
%   |_folder with masks for ROIs
%   |
%   |_folder with SPMmats
%       |_folder subject 1
%       |_folder subject 2
%       ...
% ATTENTION: In case of more than 2 masks and and more than 2 contrast
%            files: The order of masks and contrast files must fit to each
%            other. e.g. To find indiviual peaks of condition 1 in mask 1
%            and of condition 2 in mask 2 ..... !
%
% INPUT PARAMETERS:
%
% cfg
% .TR            = time of repetition
% .path2masks    = Path to the mask files
% .NrOfRuns      = numer of sessions(runs) measured per subject
% .path2files    = Path to the subjects' data files
%   and in case of several runs
%   .subdir        = cellarray with the names of the subdirectories of the
%                    data
%
% .path2SPMmat   = Path to the SPM outputfiles of the first level analyses
% .contrastname  = Cell containing the names of the contrast files used to
%                  indentify the individual peaks.
%                  ATTENTION: The order of masks and
%                  contrast files must fit to each other. e.g. indiviual
%                  peaks of condition 1 in mask 1 and of condition 2 in
%                  mask 2!
%
%
% .outputtype    = '3DAsEmbed', '3DAsTrial' or 'SingleVoxel'
%
%
% .hpfilter      = 'yes': uses a highpass filter on the fMRI data before
%                  transforming the data in the fieldtrip format.
%                  (recommended if not alreadyy done in the fMRI data
%                  preprocessing)
%   in case of cfg.hpfilter = 'yes'
%   .hpfreq      = filter frequency in Hz (eg. for strong filter you can
%                  use: (greatest differences between 2 conditions)*2 )
%
% .path4output   = path to save the output files for each subject
% .outputsuffix  = suffix for outputfiles
% .builddiff     = 'yes': build diffference from t - (t-1) for all time
%                  points to create stationary data.
% .indipeak      = 'yes': uses the individual peaks within the mask instead
%                  of the mean over the mask (default = 'no')
% .interpolate   = 'yes' or 'no' (default = 'no')
%   in case of cfg.interpolate = 'yes'
%   .interpmethod= 'liner', 'cubic' or 'spline' (default = 'spline')
%   .interpsteps = nr of interpolated data points between two real data
%                  points
% .normalize     = 'no', 'zscore'
% .selectvolumes = 'all', 'range', 'trialcut'
%   in case of 'trialcut':
%   .NrVolumes2Cut = Nr of Volumes included for each trial (including the
%                    onset)
%    .Onsets       = Cell array (nr subjects x nr of runs) containing the
%                    onsets of the condition of interest
%   in case of 'range':
%   .range       = vector including the number of first and the last volume
%                  for each run (NrOfRuns x 2)

%
% OUTPUT DATA
%
%   Data           = Fieldtrip raw data structure - containing:
%       .trial     = cell array containing the data for each trial
%       .time      = cell containing the time indices for each trial (in
%                    seconds)
%       .label     = cell containing the labels (strings) of channels
%                    included in the data
%       .fsample   = value of sampling rate (in Hertz)
%       .datatype  = 'fMRI'
%
% This Data can be used as input for TEprepare of TRENTOOL.
% The output data is not saved automatically. Save it manually if
% neccessary.
%
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation;
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY;
%
% Version 1.0 by Michael Lindner, Bonn 2011
%


% check input

if ~isfield(cfg, 'path2files');
    error('TRENTOOL ERROR: cfg.path2files has to be defined! see help')
else
    if strcmp(cfg.path2files(end),filesep)
        error('TRENTOOL ERROR: cfg.path2files must not end with a file seperator')
    end
end
if ~isfield(cfg, 'path2masks');
    error('TRENTOOL ERROR: cfg.path2masks has to be defined! see help')
else
    if strcmp(cfg.path2masks(end),filesep)
        error('TRENTOOL ERROR: cfg.path2masks must not end with a file seperator')
    end
end
if ~isfield(cfg, 'path2SPMmat');
    error('TRENTOOL ERROR: cfg.path2SPMmat has to be defined! see help')
else
    if strcmp(cfg.path2SPMmat(end),filesep)
        error('TRENTOOL ERROR: cfg.path2SMPmat must not end with a file seperator')
    end
end
if ~isfield(cfg, 'path4output');
    error('TRENTOOL ERROR: cfg.path4output has to be defined! see help')
else
    if strcmp(cfg.path4output(end),filesep)
        error('TRENTOOL ERROR: cfg.path4output must not end with a file seperator')
    end
end
if ~isfield(cfg, 'TR');
    error('TRENTOOL ERROR: cfg.TR has to be defined! see help')
end
if ~isfield(cfg, 'contrastname');
    error('TRENTOOL ERROR: cfg.contrastname has to be defined! see help')
end
if ~isfield(cfg, 'NrOfRuns');
    error('TRENTOOL ERROR: cfg.NrOfRuns has to be defined! see help')
end
if ~isfield(cfg, 'outputsuffix');
    error('TRENTOOL ERROR: cfg.outputsuffix has to be defined! see help')
end
if ~isfield(cfg, 'hpfilter');
    error('TRENTOOL ERROR: cfg.hpfilter has to be defined! see help')
end
if strcmp(cfg.hpfilter,'yes')
    if ~isfield(cfg, 'hpfreq');
        error('TRENTOOL ERROR: cfg.hpfreq has to be defined! see help')
    end
end
if ~isfield(cfg, 'outputtype');
    error('TRENTOOL ERROR: cfg.outputtype has to be defined! see help')
end

if ~isfield(cfg, 'builddiff');cfg.builddiff = 'no';end
if ~isfield(cfg, 'indipeak');cfg.indipeak = 'yes';end
if ~isfield(cfg, 'selectvolumes');
    error('TRENTOOL ERROR: cfg.selectvolumes has to be defined! see help')
end
if strcmp(cfg.selectvolumes, 'range')
    if ~isfield(cfg, 'range')
        error('TRENTOOL ERROR: cfg.range has to be defined! see help')
    end
elseif strcmp(cfg.selectvolumes, 'trialcut')
    if ~isfield(cfg, 'NrVolumes2Cut')
        error('TRENTOOL ERROR: cfg.NrVolumes2Cut has to be defined! see help')
    end
elseif strcmp(cfg.selectvolumes, 'all')
else
    error('TRENTOOL ERROR: wrong input parameter for cfg.selectvolumes! see help')
end

if ~isfield(cfg, 'normalize');cfg.normalize = 'zscore';end
if strcmp(cfg.normalize, 'no') || strcmp(cfg.normalize, 'zscore')
else
    error('TRENTOOL ERROR: wrong input for cfg.normalize! see help')
end

if ~isfield(cfg, 'interpolate'); cfg.interpolate = 'no'; end;
if strcmp(cfg.interpolate, 'yes')
    if ~isfield(cfg, 'interpmethod'); cfg.interpmethod = 'spline'; end;
    if strcmp(cfg.interpmethod, 'linear') || strcmp(cfg.interpmethod, 'cubic') || strcmp(cfg.interpmethod, 'spline')
    else
        error('TRENTOOL ERROR: wrong input for cfg.interpmethod! see help')
    end
    if ~isfield(cfg, 'interpsteps')
        error('TRENTOOL ERROR: cfg.interpsteps has to be defined! see help')
    end
    if strcmp(cfg.selectvolumes, 'trialcut')
        warning('WarnTests:convertTest',...
            '\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\nTRENTOOL WARNING: interpolation of concatenated trials could lead \nto not really existing and even maybe to wrong values in the interim of the trials\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    end
end


if strcmp(cfg.outputtype, '3DAsTrial')
    if ~isfield(cfg, 'boxsize');
        boxsize = 3;
    else
        boxsize = cfg.boxsize;
    end
elseif strcmp(cfg.outputtype, '3DAsEmbed')
    % indices of the 8 voxel connected at the edges of the middle voxel
    embindices = {[-1,-1,-1];...
        [-1,-1,1];...
        [-1,1,1];...
        [-1,1,-1];...
        [1,-1,-1];...
        [1,-1,1];...
        [1,1,1];...
        [1,1,-1]};
end



% load masks
fprintf('Load masks')

masks = dir([cfg.path2masks,'\*.nii']);
msk=cell(1,length(masks));
for ii = 1:length(masks)
    msk{ii} = load_nii([cfg.path2masks, '\' masks(ii).name]);
    if strcmp(cfg.indipeak, 'no')
        if max(size(find(msk{ii}.img==1))) ~= 1
            error(['TRENTOOL ERROR: Mask ',num2str(ii),' includes more than one voxel, see help'])
        end
    end
    cfg.roilabel{ii}=masks(ii).name;
end
fprintf(' - ok\n')

% get subject IDs
fprintf('Get subject IDs')
d = dir([cfg.path2files]);
d = d(3:end); % get rid of . and ..
NrOfSubjects=length(d);
fprintf(' - ok\n')


% load contrast images
working_directory = pwd;
if strcmp(cfg.indipeak, 'yes')
    fprintf('Load contrast images')
    if isfield(cfg, 'path2SPMmat')
        contrasts = cell(length(cfg.contrastname),NrOfSubjects);
        for subject = 1:NrOfSubjects
            for cc = 1:length(cfg.contrastname)
                cd ([cfg.path2SPMmat filesep d(subject).name filesep ])
                loadname = [cfg.contrastname{cc}, '.img'];
                contr = load_nii(loadname);
                contrasts{cc,subject} = contr;
                clear contr
            end
        end
    else
        error('TRENTOOL ERROR: cfg.path2SPMmat has to be defined!')
    end
    fprintf(' - ok\n')
end
cd(working_directory)


% Begin loop over all Subjects
for ss = 1:NrOfSubjects

    % load data
    fprintf(['Loading data: ',num2str(ss),' of ',num2str(NrOfSubjects)])
    if cfg.NrOfRuns == 1
        files = dir([cfg.path2files, filesep,d(ss).name, filesep,'*.nii']);
        for ii = 1:length(files)
            nii{ii,1} = load_nii([cfg.path2files, filesep, d(ss).name, filesep, files(ii).name]);
        end
    else
        filelength=nan(1,cfg.NrOfRuns);
        for rr = 1:cfg.NrOfRuns
            files = dir([cfg.path2files, filesep,d(ss).name, filesep,cfg.subdir{rr},'*.nii']);
            for ii = 1:length(files)
                nii{ii,rr} = load_nii([cfg.path2files, filesep,  d(ss).name, filesep, cfg.subdir{rr}, files(ii).name]);
            end
            filelength(rr)=length(files);
        end
    end
    fprintf(' - ok\n')


    % get data
    working_directory = pwd;
    cd(cfg.path2masks)

    fprintf('Reading data')
    data = cell(length(msk),cfg.NrOfRuns);
    embdata = cell(length(msk),cfg.NrOfRuns);
    if cfg.NrOfRuns == 1
        for ii = 1:length(msk)
%             dat = nan(1,length(nii));
            for jj = 1:length(nii)
                %                 if strcmp(cfg.indipeak,'yes')
                tt = contrasts{ii,ss}.img(msk{ii}.img>0);
                [xcoord,ycoord,zcoord]=ind2sub([size(contrasts{ii,ss}.img,1),size(contrasts{ii,ss}.img,2),size(contrasts{ii,ss}.img,3)],find(contrasts{ii,ss}.img==max(tt)));
                if strcmp(cfg.outputtype, '3DAsTrial')
                    voxbox = nii{jj,rr}.img(xcoord-((boxsize-1)/2):xcoord+((boxsize-1)/2),ycoord-((boxsize-1)/2):ycoord+((boxsize-1)/2),zcoord-((boxsize-1)/2):zcoord+((boxsize-1)/2));
                    voxvec = reshape(voxbox,[boxsize^3,1]);
                    dat(:,jj)=voxvec;
                elseif strcmp(cfg.outputtype, '3DAsEmbed')
                    vox = nii{jj,rr}.img(xcoord,ycoord,zcoord);
                    embvox = nan(length(embindices));
                    for extembdat =1:length(embindices)
                        embvox(extembdat) = nii{jj,rr}.img(xcoord+embindices{extembdat}(1),ycoord+embindices{extembdat}(2),zcoord+embindices{extembdat}(3));
                    end
                    dat(:,jj)=vox;
                    embdat(:,jj) = embvox;
                elseif strcmp(cfg.outputtype, 'SingleVoxel')
                    vox = nii{jj,rr}.img(xcoord,ycoord,zcoord);
                    dat(:,jj)=vox;
                end
                %                 end
            end
            data{ii,1}=dat;
            if strcmp(cfg.outputtype, '3DAsEmbed')
                embdata{ii,1}=embdat;
            end
            clear dat
        end
    else
        for ii = 1:length(msk)
            for rr = 1:cfg.NrOfRuns
                if strcmp(cfg.outputtype, '3DAsTrial')
                    dat = nan(boxsize^3,filelength(rr));
                elseif strcmp(cfg.outputtype, '3DAsEmbed')
                    dat = nan(1,filelength(rr));
                    embdat = nan(length(embindices),filelength(rr));
                elseif strcmp(cfg.outputtype, 'SingleVoxel')
                    dat = nan(1,filelength(rr));
                end
                for jj = 1:filelength(rr)
                    %                     if strcmp(cfg.indipeak,'yes')
                    tt = contrasts{ii,ss}.img(msk{ii}.img>0);
                    [xcoord,ycoord,zcoord]=ind2sub([size(contrasts{ii,ss}.img,1),size(contrasts{ii,ss}.img,2),size(contrasts{ii,ss}.img,3)],find(contrasts{ii,ss}.img==max(tt)));
                    if strcmp(cfg.outputtype, '3DAsTrial')
                        voxbox = nii{jj,rr}.img(xcoord-((boxsize-1)/2):xcoord+((boxsize-1)/2),ycoord-((boxsize-1)/2):ycoord+((boxsize-1)/2),zcoord-((boxsize-1)/2):zcoord+((boxsize-1)/2));
                        voxvec = reshape(voxbox,[boxsize^3,1]);
                        dat(:,jj)=voxvec;
                    elseif strcmp(cfg.outputtype, '3DAsEmbed')
                        vox = nii{jj,rr}.img(xcoord,ycoord,zcoord);
                        for extembdat =1:length(embindices)
                            embdat(extembdat,jj) = nii{jj,rr}.img(xcoord+embindices{extembdat}(1),ycoord+embindices{extembdat}(2),zcoord+embindices{extembdat}(3));
                        end
                        dat(:,jj)=vox;
                    elseif strcmp(cfg.outputtype, 'SingleVoxel')
                        vox = nii{jj,rr}.img(xcoord,ycoord,zcoord);
                        dat(:,jj)=vox;
                    end
                    %                     end
                end
                data{ii,rr}=dat;
                if strcmp(cfg.outputtype, '3DAsEmbed')
                    embdata{ii,rr}=embdat;
                end
                
            end
        end
    end
    clear dat
    cd(working_directory)
    fprintf(' - ok\n')

    % check nr of volumes per run of equality and cut if neccessary
    if min(filelength) ~= max(filelength)
        warning('\nTRENTOOL WARNING: Nr of data points of all runs are not identical!\nThe minimal number of data points will be used for all runs to extract the data!!\n')
    end
    for cde1 = 1:size(data,1)
        for cde2 = 1:size(data,2)
            data{cde1,cde2}=data{cde1,cde2}(:,1:min(filelength));
        end
    end
        
    
    % apply high pass filter
    if strcmp(cfg.hpfilter, 'yes')
        fprintf('High pass filter')
        for hpf1 = 1:size(data,1)
            for hpf2 = 1:size(data,2)
%                 data{hpf1,hpf2} = highpassfilter(data{hpf1,hpf2},1000/cfg.TR,cfg.hpfreq);
                data{hpf1,hpf2} = TEhighpassfilter(cfg,data{hpf1,hpf2});
                if strcmp(cfg.outputtype, '3DAsEmbed')
%                     embdata{hpf1,hpf2} = highpassfilter(embdata{hpf1,hpf2},1000/cfg.TR,cfg.hpfreq);
                    embdata{hpf1,hpf2} = TEhighpassfilter(cfg,embdata{hpf1,hpf2});
                end
            end
        end
        fprintf(' - ok\n')
    end

    
    % cut data depending on the selection type
    fprintf(['Cutting data - select volumes: ',cfg.selectvolumes])

    if strcmp(cfg.selectvolumes, 'all')
        
        % create empty matrices
        if strcmp(cfg.outputtype, '3DAsTrial')
            inrow = nan(size(data,1),boxsize^3,size(data,2));
        elseif strcmp(cfg.outputtype, '3DAsEmbed')
            inrow = nan(size(data,1),size(data,2));
            embinrow = nan(size(embdata,1),length(embindices),size(embdata,2));
        elseif strcmp(cfg.outputtype, 'SingleVoxel')
            inrow = nan(size(data,1),size(data,2));            
        end
        
        % take all data points
        for ii = 1:size(data,1)
            b = 1;
            for jj = 1:size(data,2)

                if strcmp(cfg.normalize, 'zscore')
                    if strcmp(cfg.outputtype, '3DAsTrial')
                        inrow(ii,1:boxsize^3,b:b+length(data{ii,jj})-1) = zscore(double(data{ii,jj}));
                    elseif strcmp(cfg.outputtype, '3DAsEmbed')
                        inrow(ii,b:b+length(data{ii,jj})-1) = zscore(double(data{ii,jj}));
                        embinrow(ii,1:length(embindices),b:b+length(embdata{ii,jj})-1) = zscore(double(embdata{ii,jj}));
                    elseif strcmp(cfg.outputtype, 'SingleVoxel')
                        inrow(ii,b:b+length(data{ii,jj})-1) = zscore(double(data{ii,jj}));
                    end
                    %                 elseif strcmp(cfg.normalize, 'demean')
                    %                     meandat = mean(data{ii,jj},2);
                    %                     embmeandat = mean(embdata{ii,jj},2);
                    %                     for kkk = 1:boxsize^3
                    %                         inrow(ii,kkk,b:b+length(data{ii,jj})-1) = data{ii,jj}(kkk,:) - meandat(kkk);
                    %                         embinrow(ii,kkk,b:b+length(embdata{ii,jj})-1) = embdata{ii,jj}(kkk,:) - embmeandat(kkk);
                    %                     end
                else
                    if strcmp(cfg.outputtype, '3DAsTrial')
                        inrow(ii,1:boxsize^3,b:b+length(data{ii,jj})-1) = data{ii,jj};
                    elseif strcmp(cfg.outputtype, '3DAsEmbed')
                        inrow(ii,b:b+length(data{ii,jj})-1) = data{ii,jj};
                        embinrow(ii,1:length(embindices),b:b+length(embdata{ii,jj})-1) = embdata{ii,jj};
                    elseif strcmp(cfg.outputtype, 'SingleVoxel')
                        inrow(ii,b:b+length(data{ii,jj})-1) = data{ii,jj};
                    end
                end
                b=b+length(data{ii,jj});
            end
        end

    elseif strcmp(cfg.selectvolumes, 'range')
        
        % create empty matrices
        if strcmp(cfg.outputtype, '3DAsTrial')
            inrow = nan(size(data,1),boxsize^3,(cfg.range(1,2)-cfg.range(1,1))+ (cfg.range(2,2)-cfg.range(2,1)) +2);
        elseif strcmp(cfg.outputtype, '3DAsEmbed')
            inrow = nan(size(data,1),(cfg.range(1,2)-cfg.range(1,1))+ (cfg.range(2,2)-cfg.range(2,1)) +2);
            embinrow = nan(size(embdata,1),length(embindices),(cfg.range(1,2)-cfg.range(1,1))+ (cfg.range(2,2)-cfg.range(2,1)) +2);
        elseif strcmp(cfg.outputtype, 'SingleVoxel')
            inrow = nan(size(data,1),(cfg.range(1,2)-cfg.range(1,1))+ (cfg.range(2,2)-cfg.range(2,1)) +2);
        end
        
        % cut range of data
        for ii = 1:size(data,1)
            b = 1;
            for jj = 1:size(data,2)

                if strcmp(cfg.normalize, 'zscore')
                    if strcmp(cfg.outputtype, '3DAsTrial')
                        inrow(ii,1:boxsize^3,b:b+length(cfg.range(jj,1):cfg.range(jj,2))-1) = zscore(double(data{ii,jj}(cfg.range(jj,1):cfg.range(jj,2))));
                    elseif strcmp(cfg.outputtype, '3DAsEmbed')
                        inrow(ii,b:b+length(cfg.range(jj,1):cfg.range(jj,2))-1) = zscore(double(data{ii,jj}(cfg.range(jj,1):cfg.range(jj,2))));
                        embinrow(ii,1:length(embindices),b:b+length(cfg.range(jj,1):cfg.range(jj,2))-1) = zscore(double(embdata{ii,jj}(cfg.range(jj,1):cfg.range(jj,2))));
                    elseif strcmp(cfg.outputtype, 'SingleVoxel')
                        inrow(ii,b:b+length(cfg.range(jj,1):cfg.range(jj,2))-1) = zscore(double(data{ii,jj}(cfg.range(jj,1):cfg.range(jj,2))));
                    end
                    %                 elseif strcmp(cfg.normalize, 'demean')
                    %                     inrow(ii,1:boxsize^3,b:b+length(cfg.range(jj,1):cfg.range(jj,2))-1) = data{ii,jj} - mean(data{ii,jj}(cfg.range(jj,1):cfg.range(jj,2)),2);
                    %                     embinrow(ii,1:boxsize^3,b:b+length(cfg.range(jj,1):cfg.range(jj,2))-1) = embdata{ii,jj} - mean(embdata{ii,jj}(cfg.range(jj,1):cfg.range(jj,2)),2);
                else
                    if strcmp(cfg.outputtype, '3DAsTrial')
                        inrow(ii,1:boxsize^3,b:b+length(cfg.range(jj,1):cfg.range(jj,2))-1) = data{ii,jj}(cfg.range(jj,1):cfg.range(jj,2));
                    elseif strcmp(cfg.outputtype, '3DAsEmbed')
                        inrow(ii,b:b+length(cfg.range(jj,1):cfg.range(jj,2))-1) = data{ii,jj}(cfg.range(jj,1):cfg.range(jj,2));
                        embinrow(ii,1:length(embindices),b:b+length(cfg.range(jj,1):cfg.range(jj,2))-1) = embdata{ii,jj}(cfg.range(jj,1):cfg.range(jj,2));
                    elseif strcmp(cfg.outputtype, 'SingleVoxel')
                        inrow(ii,b:b+length(cfg.range(jj,1):cfg.range(jj,2))-1) = data{ii,jj}(cfg.range(jj,1):cfg.range(jj,2));
                    end
                end
                b=b+length(cfg.range(jj,1):cfg.range(jj,2));
            end
        end

    elseif strcmp(cfg.selectvolumes, 'trialcut')

        % check nr of onsets for equality
        NOm = nan(size(cfg.Onsets,1),size(cfg.Onsets,2));
        for s1 = 1:size(cfg.Onsets,1);
            for s2 = 1:size(cfg.Onsets,2);
                NOm(s1,s2) = length(cfg.Onsets{s1,s2});
            end
        end
        clear s1 s2
        NO=sum(NOm,2);
        if min(NO) ~= max(NO)
            error('TRENTOOL ERROR: Number of Onsets must be equal for all')
        end

        % calculate onsets and offset volumes
        fprintf(' - (Calculating on- and offsets')
        Onsets = cell(1,size(cfg.Onsets,2));
        Offsets = cell(1,size(cfg.Onsets,2));
        for ii = 1:size(cfg.Onsets,2)
            Onsets{ii} = floor(cfg.Onsets{ss,ii}/(cfg.TR/1000));
            Onsets{ii}(find(Onsets{ii}==0))=1;
            Offsets{ii} = Onsets{ii}+cfg.NrVolumes2Cut-1;
        end
        fprintf(' - ok ) ')
        
        %create empty matrices
        if strcmp(cfg.outputtype, '3DAsTrial')
            inrow = nan(size(data,1), boxsize^3, min(NO)*cfg.NrVolumes2Cut);
        elseif strcmp(cfg.outputtype, '3DAsEmbed')
            inrow = nan(size(data,1), min(NO)*cfg.NrVolumes2Cut);
            embinrow = nan(size(data,1), length(embindices), min(NO)*cfg.NrVolumes2Cut);
        elseif strcmp(cfg.outputtype, 'SingleVoxel')
            inrow = nan(size(data,1), min(NO)*cfg.NrVolumes2Cut);
        end
        
        % Cut trial data 
        for kk = 1:size(data,1)
            b = 1;
            for oo = 1:length(Onsets)
                a=1;
                for ff = 1:length(Onsets{oo})
                    if strcmp(cfg.normalize, 'zscore')
                        if strcmp(cfg.outputtype, '3DAsTrial')
                            vec(1:boxsize^3,a:a+cfg.NrVolumes2Cut-1) = zscore(double(data{kk,oo}(:,Onsets{oo}(ff):Offsets{oo}(ff))));
                        elseif strcmp(cfg.outputtype, '3DAsEmbed')
                            vec(a:a+cfg.NrVolumes2Cut-1) = zscore(double(data{kk,oo}(:,Onsets{oo}(ff):Offsets{oo}(ff))));
                            embvec(1:length(embindices),a:a+cfg.NrVolumes2Cut-1) = zscore(double(embdata{kk,oo}(:,Onsets{oo}(ff):Offsets{oo}(ff))));
                        elseif strcmp(cfg.outputtype, 'SingleVoxel')
                            vec(a:a+cfg.NrVolumes2Cut-1) = zscore(double(data{kk,oo}(:,Onsets{oo}(ff):Offsets{oo}(ff))));
                        end
                        %                     elseif strcmp(cfg.normalize, 'demean')
                        %                         vec(1:boxsize^3,a:a+cfg.NrVolumes2Cut-1) = data{ii,oo}(:,Onsets{oo}(ff):Offsets{oo}(ff)) - mean(data{ii,oo}(Onsets{oo}(ff):Offsets{oo}(ff)));
                    else
                        if strcmp(cfg.outputtype, '3DAsTrial')
                            vec(1:boxsize^3,a:a+cfg.NrVolumes2Cut-1) = data{kk,oo}(:,Onsets{oo}(ff):Offsets{oo}(ff));
                        elseif strcmp(cfg.outputtype, '3DAsEmbed')
                            vec(a:a+cfg.NrVolumes2Cut-1) = data{kk,oo}(:,Onsets{oo}(ff):Offsets{oo}(ff));
                            embvec(1:length(embindices),a:a+cfg.NrVolumes2Cut-1) = embdata{kk,oo}(:,Onsets{oo}(ff):Offsets{oo}(ff));
                        elseif strcmp(cfg.outputtype, 'SingleVoxel')
                            vec(a:a+cfg.NrVolumes2Cut-1) = data{kk,oo}(:,Onsets{oo}(ff):Offsets{oo}(ff));
                        end
                    end
                    a=a+cfg.NrVolumes2Cut;
                end
                if strcmp(cfg.outputtype, '3DAsTrial')
                    inrow(kk,1:boxsize^3,b:b+size(vec,2)-1) = vec;
                elseif strcmp(cfg.outputtype, '3DAsEmbed')
                    inrow(kk,b:b+size(vec,2)-1) = vec;
                    embinrow(kk,1:length(embindices),b:b+size(vec,2)-1) = embvec;
                elseif strcmp(cfg.outputtype, 'SingleVoxel')
                    inrow(kk,b:b+size(vec,2)-1) = vec;
                end
                b=b+size(vec,2);
                clear vec embvec
            end
        end
    end
    clear a b oo ii
    fprintf(' - ok\n')

    % interpolate data
    if strcmp(cfg.interpolate, 'yes')
        fprintf('interpolate data')
        
        if strcmp(cfg.outputtype, '3DAsTrial')
            x = 1:size(inrow,3);
            xx = 1:1/(cfg.interpsteps):size(inrow,3);
            inrow_interp=nan(size(inrow,1),boxsize^3,size(xx,2));
        elseif strcmp(cfg.outputtype, '3DAsEmbed')
            x = 1:size(inrow,2);
            xx = 1:1/(cfg.interpsteps):size(inrow,2);
            inrow_interp=nan(size(inrow,1),size(xx,2));
            embinrow_interp=nan(size(embinrow,1),length(embindices),size(xx,2));
        elseif strcmp(cfg.outputtype, 'SingleVoxel')
            x = 1:size(inrow,2);
            xx = 1:1/(cfg.interpsteps):size(inrow,2);
            inrow_interp=nan(size(inrow,1),size(xx,2));
        end
        for ii = 1:size(inrow,1)
            if strcmp(cfg.outputtype, '3DAsTrial')
                for mm = 1:boxsize^3
                    inrow_interp(ii,mm,:) = interp1(x,squeeze(inrow(ii,mm,:))',xx,cfg.interpmethod);
                end
            elseif strcmp(cfg.outputtype, '3DAsEmbed')
                inrow_interp(ii,:) = interp1(x',squeeze(inrow(ii,:))',xx,cfg.interpmethod);
                for mm = 1:length(embindices)
                    embinrow_interp(ii,mm,:) = interp1(x,squeeze(embinrow(ii,mm,:))',xx,cfg.interpmethod);
                end
            elseif strcmp(cfg.outputtype, 'SingleVoxel')
                inrow_interp(ii,:) = interp1(x',squeeze(inrow(ii,:))',xx,cfg.interpmethod);
            end
        end
        clear inrow embinrow
        inrow=inrow_interp;
        if strcmp(cfg.outputtype, '3DAsEmbed')
            embinrow=embinrow_interp;
        end
        clear inrow_interp embinrow_interp x x ii mm
        fprintf(' - ok\n')
    end

    % build difference
    if strcmp(cfg.builddiff, 'yes')
        fprintf('Build differences')
        if strcmp(cfg.outputtype, '3DAsTrial')
            inrow = inrow(:,:,2:end)-inrow(:,:,1:end-1);
        elseif strcmp(cfg.outputtype, '3DAsEmbed')
            inrow = inrow(:,2:end)-inrow(:,1:end-1);
            embinrow = embinrow(:,:,2:end)-embinrow(:,:,1:end-1);
        elseif strcmp(cfg.outputtype, 'SingleVoxel')
            inrow = inrow(:,2:end)-inrow(:,1:end-1);
        end
        fprintf(' - ok\n')
    end
    
    % test stationarity??
    
    %build trials
    fprintf('Build trials')
    if strcmp(cfg.outputtype, '3DAsTrial')
        for vt = 1:boxsize^3
            Data.trial{vt}=squeeze(inrow(:,vt,:));
        end
    elseif strcmp(cfg.outputtype, '3DAsEmbed')
        Data.trial{ss}=inrow;
        Data.Data4Embedding{ss}=embinrow;
    elseif strcmp(cfg.outputtype, 'SingleVoxel')
        Data.trial{ss}=inrow;
    end
    fprintf(' - ok\n')

    % prepare out put if '3DAsTrial'
    if strcmp(cfg.outputtype, '3DAsTrial')
        % prepare output
        fprintf('Prepare output Data')
        for tttt=1:boxsize^3
            Data.time{tttt}=0:cfg.TR/1000:(size(inrow,3)-1)*cfg.TR/1000;

        end
        Data.label=cfg.roilabel;
        Data.fsample = 1000/cfg.TR;
        Data.datatype = 'fMRI';
        Data.outputtype = cfg.outputtype;
        fprintf(' - ok\n')

        % save file
        try
            cd(cfg.path4output)
        catch
            mkdir(cfg.path4output)
            cd(cfg.path4output)
        end
        savename = [cfg.path4output,filesep,d(ss,1).name,cfg.outputsuffix];
        fprintf(['save data file: ',d(ss,1).name,cfg.outputsuffix])
        save(savename,'Data');
        fprintf(' - ok\n')
        clear Data savename
    end
    if ss == 1
        s4t = size(inrow,2);
    end
    clear data dat Onsets Offsets vec inrow embinrow nii 
end

if strcmp(cfg.outputtype, '3DAsEmbed') || strcmp(cfg.outputtype, 'SingleVoxel') 
        try
            cd(cfg.path4output)
        catch
            mkdir(cfg.path4output)
            cd(cfg.path4output)
        end
    fprintf('Prepare output Data')
    for tttt=1:NrOfSubjects
        Data.time{tttt}=0:cfg.TR/1000:(s4t-1)*cfg.TR/1000;

    end
    Data.label=cfg.roilabel;
    Data.fsample = 1000/cfg.TR;
    Data.datatype = 'fMRI';
    Data.outputtype = cfg.outputtype;
    fprintf(' - ok\n')
    
    % save file
    if strcmp(cfg.outputtype, '3DAsEmbed')
        savename = [cfg.path4output,filesep,'SubjAsTrials_3DAsEmbed',cfg.outputsuffix];
    elseif strcmp(cfg.outputtype, 'SingleVoxel')
        savename = [cfg.path4output,filesep,'SubjAsTrials_SingleVoxel',cfg.outputsuffix];
    end
    fprintf(['save data file: ',savename])
    save(savename,'Data');
    fprintf(' - ok\n')

    fprintf('Thanks for using this function!\n\ndone\n')
end
