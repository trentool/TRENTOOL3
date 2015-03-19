function TRENTOOL2BrainNet(cfg,TEpermtest)

% FUNCTION TRENTOOL2BrainNet(cfg,TEpermtest,MNIcoord,labels,filename,nodeCol,nodeSize,edgeCol)
%
% The function converts TRENTOOL outputs to *.node and *.edge files
% readable by BrainNet Viewer (Xia, 2013, PlosONE, 
% http://www.nitrc.org/projects/bnv/).
%
% * INPUT
%   TEpermtest = structure containing TE results, returned by TRENTOOL
%                functions TEsurrogatestats.m, TEsurrogatestats_ensemble.m 
%                or InteractionDelayReconstruction_analyze.m
%
%  cfg. 	 Configuration structure with fields
%   MNIcoord   = MNI coordinates (x,y,z) of the sources, array with size 
%                [N 3], where N is the number of sources
%   labels     = source labels, cell array with size [N 1] - Not that
%                sources shouldn't contain spaces (causes BrainNet to
%                crash), the function will remove any spaces from the 
%                labels cell array. WARNING: You may use the labels
%                provided in TEpermtest.cfg.channel, if you use different
%                labels, make sure, they are in the same order as the
%                channels in TEpermtest.cfg.channel!
%   filename   = filename for output files *.node and *.edge, you may
%                specify a filename only (string), which causes the 
%                function to save both files to the current folder, or you 
%                can specify a whole filepath ([filepath]/[filename])
%
%
%   THE FOLLOWING INPUTS ARE OPTIONAL:
%
%  cfg.
%   nodeCol    = BrainNet gives you the option to color nodes according to
%                the values in this vector (size [N 1]), see BrainNet
%                Manual. If you want the nodes to have the same color,
%                provide no vector or a vector containing ones (default).
%   nodeSize   = BrainNet gives you the option to size nodes according to
%                the values in this vector (size [N 1]), see BrainNet
%                Manual. If you want the nodes to have the same size,
%                provide no vector or a vector containing ones (default).
%   edgeCol    = BrainNet gives you the option to color edges according to
%                the values in this vector (size [N 1]), see BrainNet
%                Manual. If you want the edges to have the same color,
%                provide no vector or a vector containing ones (default).
%
%
% PW, Frankfurt, 2013

% Changelog
% 30/10/2014: PW fixed a bug in the edge coloring


%% check input and get parameters
MNIcoord = cfg.MNIcoord;
labels   = cfg.labels;
filename = cfg.filename;

N = size(labels,1);
E = size(TEpermtest.sgncmb,1);
signValue = 2; % use uncorrected TE value

if nargin < 2
    error('Not enough input arguments!')
end
   
if ~isfield(cfg,'nodeCol');  nodeCol  = ones(N,1); end;
if ~isfield(cfg,'nodeSize'); nodeSize = ones(N,1); end;
if ~isfield(cfg,'edgeCol');  
    edgeCol  = ones(E,1); 
else
    edgeCol  = cfg.edgeCol;
end;

% check, if there are blanks in the labels and delte them, otherwise 
% BrainNet chrashes
for i=1:N;
    % find blanks in labels
    blanks = strfind(labels{i},' ');
    
    % if there are blanks, delete them
    if ~isempty(blanks)
        for j=1:length(blanks)
            labels{i}(blanks(j)) = [];
        end
    end
end

%% write node file
fid = fopen([filename '.node'],'w');

for i = 1:size(MNIcoord,1)
    fprintf(fid,'%.0f\t%.0f\t%.0f\t',MNIcoord(i,:));
    fprintf(fid,'%.0f\t%.0f\t',nodeCol(i),nodeSize(i));
    fprintf(fid,'%s\n',labels{i});
end

fclose(fid);

%% write edge file

adjacencyMatrix = zeros(N);

for i=1:E
    
    % check if TE is significant for signalcombination i
    if TEpermtest.TEpermvalues(i,signValue)
    
        % find the indices of the channels in signalcombination i
        ind1 = find(strcmp(TEpermtest.sgncmb{i,1},TEpermtest.cfg.channel));
        ind2 = find(strcmp(TEpermtest.sgncmb{i,2},TEpermtest.cfg.channel));
    
        adjacencyMatrix(ind1,ind2) = edgeCol(i);
    end
end


fid = fopen([filename '.edge'],'w');

for i=1:N
    for j=1:N
        fprintf(fid,'%.0f\t',adjacencyMatrix(i,j));
    end
    fprintf(fid,'\n');
end

fclose(fid);
disp(['done -  files have been saved to ''' filename '.edge'' and ''' filename '.node''.']);