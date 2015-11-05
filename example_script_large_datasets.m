
%% set paths
addpath('~/repos/TRENTOOL3/') % note: this script assumes, you have already compiled the GPU mex-files using 'install.m'
addpath('~/TRENTOOL_gpu_functions/')
addpath('~/fieldtrip-20150928/');
ft_defaults;

%%
data_path = '~/repos/TRENTOOL3_exampledata/Mooney/';
output_path = '~/TRENTOOL_dev/test_current_version/results/';
files = dir(sprintf('%sMooney_ex*', data_path));
load(sprintf('%s%s', data_path, files(1).name))

%% define cfg for TEprepare.m

cfgTEP = [];

% data
cfgTEP.toi                 = [min(VChannelDataOut.time{1,1}),max(VChannelDataOut.time{1,1})]; % time of interest
cfgTEP.sgncmb              = {... % channels to be analyzed
    VChannelDataOut.label{1} VChannelDataOut.label{2}; ...  
    VChannelDataOut.label{1} VChannelDataOut.label{3}; ...
    VChannelDataOut.label{1} VChannelDataOut.label{4}};

% scanning of interaction delays u
cfgTEP.predicttimemin_u    = 3;		  % minimum u to be scanned
cfgTEP.predicttimemax_u    = 7;	  % maximum u to be scanned
cfgTEP.predicttimestepsize = 2; 	  % time steps between u's to be scanned

% use ensemble method
cfgTEP.ensemblemethod = 'yes';

% ACT estimation and constraints on allowed ACT(autocorelation time)
cfgTEP.actthrvalue = 40;   % threshold for ACT
cfgTEP.maxlag      = 100;
cfgTEP.minnrtrials = 15;   % minimum acceptable number of trials

% optimizing embedding
cfgTEP.optimizemethod ='ragwitz';  % criterion used
cfgTEP.ragdim         = 2:8;       % criterion dimension
cfgTEP.ragtaurange    = [0.2 0.4]; % range for tau
cfgTEP.ragtausteps    = 15;        % steps for ragwitz tau steps
cfgTEP.repPred        = 100;       % size(data.trial{1,1},2)*(3/4);

% kernel-based TE estimation
cfgTEP.flagNei = 'Mass' ;           % neigbour analyse type
cfgTEP.sizeNei = 4;                 % neigbours to analyse


%% define cfg for TEsurrogatestats_ensemble.m

cfgTESS = [];

% use individual dimensions for embedding
cfgTESS.optdimusage = 'indivdim';
cfgTESS.surrogatetype  = 'trialperm';   % surrogate type
cfgTESS.numpermutation = 1;             % no. surrogate data sets

% GPU specifications
cfgTESS.GPUmemsize     = 4200;
cfgTESS.numthreads     = 512;
cfgTESS.maxgriddim     = 65535;

cfgTESS.extracond      = 'Faes_Method'; % control for volume conduction
cfgTESS.shifttest      = 'no';
cfgTESS.MIcalc = 0; % don't calculate MI additionally to TE
cfgTESS.fileidout  = sprintf('%sLorenzdata_1->2_ensemble_', output_path); % results file name


%%

for subj = 1:length(files)    
    load(sprintf('%s%s', data_path, files(subj).name))
    InteractionDelayReconstruction_calculate(cfgTEP,cfgTESS,VChannelDataOut);   
end

%% 
output_files = dir(sprintf('%sMooney_ex_*_RAG4_TGA_opt_u*_TEpermtest_output.mat', output_path));

run('~/repos/TRENTOOL3/private/TEsetRandStream.m')
numpart = round(47+rand(1,10)*10);
charpart = char(round(65+rand(1,10)*26));
code = [numpart(1:5),charpart(1:5),numpart(6:10),charpart(6:10)];

% prepare output structure
groupprepare = [];
groupprepare.code             = code;
groupprepare.predicttimevec_u = ...
    cfgTEP.predicttimemin_u:cfgTEP.predicttimestepsize:cfgTEP.predicttimemax_u;

for subj = 1:length(output_files)
    load(sprintf('%s%s', output_path, output_files(subj).name))   
    TEpermtest.groupprepare = groupprepare;
    TEpermtest_sur = TEpermtest;
    TEpermtest_sur.TEmat = TEpermtest_sur.TEmat_sur;
    save(sprintf('%s%s', output_path, output_files(subj).name), ...
        'TEpermtest');
    clear TEpermtest;
    TEpermtest = TEpermtest_sur;
    save(sprintf('%s%s_surrogate.mat', output_path, output_files(subj).name(1:end-4)), ...
        'TEpermtest');
end

%%

% collect data in a cell array of file names
cd(output_path)
orig_files = dir( ...
    sprintf('%sMooney_ex_*_RAG4_TGA_opt_u*_TEpermtest_output.mat', output_path));
sur_files  = dir( ...
    sprintf('%sMooney_ex_*_RAG4_TGA_opt_u*_TEpermtest_output_surrogate.mat', output_path));
stats_files = {orig_files(:).name sur_files(:).name};

% define analysis parameters
cfgGSTAT = [];
cfgGSTAT.design    = [1 2 3 1 2 3; 1 1 1 2 2 2];
cfgGSTAT.uvar      = 1;
cfgGSTAT.ivar      = 2;
cfgGSTAT.fileidout      = 'Mooney_new_stats';
cfgGSTAT.permstatstype  = 'depsamplesT';

% run group statistics
TEgroup_stats(cfgGSTAT, stats_files);
