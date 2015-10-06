%% set paths

addpath('TRENTOOL3')
addpath('fieldtrip-20150928');
ft_defaults;

%% define data paths

OutputDataPath = '~/results/';
InputDataPath  = '~/TRENTOOL3_exampledata/Lorenz_2_systems/lorenz_1-2_45ms.mat';

load(InputDataPath);


%% define cfg for TEprepare.m

cfgTEP = [];


% data
cfgTEP.toi                 = [min(data.time{1,1}),max(data.time{1,1})]; % time of interest
cfgTEP.sgncmb              = {'A1' 'A2'};  % channels to be analyzed

% scanning of interaction delays u
cfgTEP.predicttimemin_u    = 40;      % minimum u to be scanned
cfgTEP.predicttimemax_u    = 50;	  % maximum u to be scanned
cfgTEP.predicttimestepsize = 1; 	  % time steps between u's to be scanned

% estimator
cfgTEP.TEcalctype  = 'VW_ds'; % use the new TE estimator (Wibral, 2013)

% ACT estimation and constraints on allowed ACT(autocorelation time)
cfgTEP.actthrvalue = 100;   % threshold for ACT
cfgTEP.maxlag      = 1000;
cfgTEP.minnrtrials = 15;   % minimum acceptable number of trials

% optimizing embedding
cfgTEP.optimizemethod ='ragwitz';  % criterion used
cfgTEP.ragdim         = 2:9;       % criterion dimension
cfgTEP.ragtaurange    = [0.2 0.4]; % range for tau
cfgTEP.ragtausteps    = 5;        % steps for ragwitz tau steps
cfgTEP.repPred        = 100;      % size(data.trial{1,1},2)*(3/4);

% kernel-based TE estimation
cfgTEP.flagNei = 'Mass' ;           % neigbour analyse type
cfgTEP.sizeNei = 4;                 % neigbours to analyse


%% define cfg for TEsurrogatestats_ensemble.m

cfgTESS = [];

% use individual dimensions for embedding
cfgTESS.optdimusage = 'indivdim';

% statistical and shift testing
cfgTESS.tail           = 1;
cfgTESS.numpermutation = 5e4;
cfgTESS.shifttesttype  ='TEshift>TE';
cfgTESS.surrogatetype  = 'trialshuffling';

% results file name
cfgTESS.fileidout  = strcat(OutputDataPath,'Lorenzdata_1->2_');

%% calculation - scan over specified values for u

TGA_results = InteractionDelayReconstruction_calculate(cfgTEP,cfgTESS,data);

save([OutputDataPath 'Lorenz_1->2_TGA_results.mat'],'TGA_results');


%% optional: perform a post hoc correction for cascade effects and simple common drive effects

cfgGA = [];

cfgGA.threshold = 3;
cfgGA.cmc       = 1;

TGA_results_GA = TEgraphanalysis(cfgGA,TGA_results);

save([OutputDataPath 'Lorenz_1->2_TGA_results_analyzed_GA.mat'],'TGA_results_GA');



%% plotting

load('exampledata/Lorenz_layout.mat');

cfgPLOT = [];

cfgPLOT.layout        = lay_Lorenz; 		% see fieldtrip's ft_prepare_layout.m
cfgPLOT.electrodes    = 'highlights';
cfgPLOT.statstype     = 1;   		% 1: corrected; 2:uncorrected; 3: 1-pval; 4:rawdistance
cfgPLOT.alpha         = 0.05;
cfgPLOT.arrowpos      = 1;
cfgPLOT.showlabels    = 'yes';
cfgPLOT.electrodes    = 'on';
cfgPLOT.hlmarker      = 'o';
cfgPLOT.hlcolor       = [0 0 0];
cfgPLOT.hlmarkersize  = 4;
cfgPLOT.arrowcolorpos = [1 0 0];

figure; 
TEplot2D(cfgPLOT,TGA_results_GA)



