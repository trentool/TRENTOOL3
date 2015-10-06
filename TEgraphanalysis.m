function [data_paths,cfg] = TEgraphanalysis(cfg,data)

% TEGRAPHANALYSIS: Detects potentially spurious edges in a graph
% constructed from single subject or single condition TE data. Note that
% the delay times u have to be integer values.
%
% TEgraphanalysis is a wrapper-function, calls TEdfs, TEdyn and 
% TEbacktracking to detect alternative paths for any edge in the TE-graph.
%
% The function constructs a graph from the provided data, such that
%    - edgess are defined by data.sgncmb (significant interactions only), 
%    - edge-weights are defined by delay times from data.TEpermvalues,
%    - vertices are enumerated according to their appearance in data.sgncmb
%    - at this point the function considers interactions SIGNIFICANT AT THE 
%      PRESCRIBED ALPHA LEVEL only (no correction for multiple comparison)!
%
% Than the function iteratively
%   - removes an edge from the graph (the weight of this edge is defined as 
%     w_crit)
%   - looks for alternative paths for this edge by running TEdyn
%   - if an alternative path exists, it is reconstructed by TEbacktracking
%     (note that TEbacktracking aborts after a certain number of paths 
%     is found to avoid excessive running times, in this case 'Too many 
%     paths. Return.' is displayed in the console, paths that have been
%     found up to that point are kept and the current edge is flagged)
%
% Alternative paths are collected for all edges. Finally, spurious edges
% are flagged by calling TEflagedges.
% See the reference information for a more detailed description.
%
%
% * REFERENCE INFORMATION
%
%   - graph algorithm
%         - Bsc Thesis Patricia Wollstadt 
%           (email: p.wollstadt@stud.uni-frankfurt.de)
%
%
% * DEPENDENCIES
%     - The functions
%         - TEdfs
%         - TEdyn
%         - TEbacktracking
%         - TEconsoleoutput
%     - FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip, Copyright 
%       (C) 2004-2008, Robert Oostenveld (GNU General Public License)
%         - ft_progress
%
%
% * INPUT PARAMETERS
%
%   cfg: The configuration MUST contain:
%
%       .threshold    = scalar value (in ms): tolerance that is used to 
%                       define the reconstruction interval around w_crit - 
%                       all paths that have a summed weight that falls 
%                       within this interval are considered an alternative 
%                       path
%       .cmc          = tells the function whether to use links that are
%                       significant after correction for multiple
%                       comparison (1) or links that are significant at the
%                       original alpha level (0)
%       .verbosity    = set the verbosity of console output (see 'help
%                       TEconsoleoutput', default: 'info_minor')
%
%   data
%       .TEprepare.channelcombi = 2xN matrix that defines analyzed channel
%                                 combinations by indices, indexing the
%                                 field channellabel
%       .TEprepare.channellabel = list of strings holding channel labels
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
%                           6 - delay times u
%
%
% * OUTPUT PARAMETERS 
%   data
%       .TEpermvalues = matrix with size channelpair x 6 (for the exact
%                       specification see INPUT PARAMETERS, if an 
%                       alternative path was found the following changes 
%                       are made for the respectice channelpair:
%                           1 - p-value is set to 1
%                           2 - significance at the prescribed alpha level 
%                               is set to 0
%                           3 - significance after correction for 
%                               multiple comparison is set to 0
%                           4 - mean difference is set to NaN
%                           5 - is set to 2/3/4 according to the type of
%                               spurious interaction:   
%                               2 = cascade effect
%                               3 = cascade effect triangle
%                               4 = common drive link triangle
%                           6 - delay times are set to 0
%
%       .graphanalysis   = contains information on the constructed graph as
%                          n_vertices = number of vertices, 
%                          n_edges    = number of edges
%                          density    = graph density, defined as
%                                       dens = E/(V*(V-1))
%                                       V = n_vertices and E = n_edges
%                          cfg        = user provided analysis parameters
%                                       (see INPUT PARAMETERS)
%
% PW - 07/09/2012
%
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation;
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY;
%
% Version 1.21 by Patricia Wollstadt
% Frankfurt 2012

% CHANGELOG
%
% 2012-04-27: PW added exception for graphs with only two nodes,
% computation is now aborted with a hint, that this doesn't make sense
%
% 2012-06-06: PW changed naming of vertices -> individual sources are
% enumerated, enumeration is used for the naming of vertices
%
% 2012-06-08: PW added documentaion
%
% 2012-07-04: PW added documentaion and references and changed command-line 
%             feedback (uses now ft_progress)
%
% 2012-07-04: NP minor changes
%
% 2013-02-13: PW added documentation for the case of too many paths 
%             during backtracking
%
% 2013-03-15: PW fixed an error in the construction of the predecessor
%             structure/adjacency matrix construction (crashed if a node
%             had no predecessors)
%
% 2014-03-21: PW added the list of triangles to the output
%
% 2015-03-20: PW changed the output structure if no paths are found

%% define logging levels
LOG_INFO_MAJOR   = 1;
LOG_INFO_MINOR   = 2;
LOG_DEBUG_COARSE = 3;
LOG_DEBUG_FINE   = 4;

% check if a threshold is provided
if ~isfield(data.TEprepare.cfg,'verbosity')
    cfg.verbosity = 'info_minor'; 
else
    cfg.verbosity = data.TEprepare.cfg.verbosity;
end;

msg = '################### CORRECTING FOR POTENTIALLY SPURIOUS EDGES';
TEconsoleoutput(cfg.verbosity, msg, LOG_INFO_MAJOR);


% check if a threshold is provided
if isfield(cfg,'threshold')
	threshold = cfg.threshold;
else
	error('No threshold defined');
end;

% check if fieldtrip version is new enough
if exist('ft_progress','file')==0;
    error('You have no current fieldtrip version in your path, that provides the function ft_progress. Please update to a version fieldtrip-201201xx or higher.')
end;

% check if edge weights are integer numbers
if sum(mod(data.TEpermvalues(:,6),1))>0, error('Delay times have to be integer values!'), end;  

% find edges and weights from input data, decide whether to use links
% significant after correction for multiple comparison

if ~isfield(cfg,'cmc')
    error('Please specify the use of corrected or uncorrected significance in ''cfg.cmc''! See help.');
end
if cfg.cmc == 1
    link_ind = data.TEpermvalues(:,3) == 1;    
    msg = 'Using links that are significant after correction for multiple comparisons';
    TEconsoleoutput(cfg.verbosity, msg, LOG_INFO_MINOR);
elseif cfg.cmc == 0
    link_ind = data.TEpermvalues(:,2) == 1;
    msg = 'Using links that are significant without correction for multiple comparisons';
    TEconsoleoutput(cfg.verbosity, msg, LOG_INFO_MINOR);
end
weights        = data.TEpermvalues(link_ind,6);
edges4analysis = data.TEprepare.channelcombi(link_ind,:);       % use numeric representation of channel combis
edges_original = data.TEprepare.channelcombi;                   %   rather than labels (strings)

% remember graph info
labels_vertices = data.TEprepare.channellabel;
n_vertices      = length(labels_vertices);
n_edges         = length(edges4analysis);


% generate output structure, graph-related info goes into a seperate substructure
msg = sprintf('no of edges: %d , no of vertices: %d', n_edges, n_vertices);
TEconsoleoutput(cfg.verbosity, msg, LOG_INFO_MINOR);
graphanalysis = [];
graphanalysis.n_edges = n_edges;
graphanalysis.n_vertices = n_vertices;
graphanalysis.density  = getDensity(n_vertices, n_edges);
graphanalysis.TEpermvalues_old = data.TEpermvalues;

% collects all alternative paths within the reconstruction interval 
% for all edges of the graph, contains: 
% [edge number] [starting node] [target node] [number of altpaths] [TEbacktracking output]
all_paths = {};


%% check if graph is big enough for graph analysis, else return

if n_vertices < 3 || n_edges < 3;
    warning(['\nTRENTOOL WARNING: The input graph has %d nodes and %d ' ...
        'edges. Graphanalysis does not work for graphs with less ' ... 
        'than 3 nodes or less than 3 edges! Return...'], ...
        n_vertices, n_edges);
        
    % add graph info to datastructure
    data_paths                   = data;
    data_paths.graphanalysis     = graphanalysis;
    data_paths.graphanalysis.cfg = cfg;
    
    % update TEsteps
    if ~isfield(data,'TEsteps')
        data_paths.TEsteps = 'GA';
    else
        data_paths.TEsteps = strcat(data.TEsteps,'_GA');
    end
    return;
end;


%% find alternative paths for all neighbours
% init progress bar
if ~strcmp(cfg.verbosity, 'none')
    fprintf('\n')
    stack = dbstack;
    msg = [ ...
            repmat('   ', 1, length(stack)-1) ...
            stack(1).file ...
            ' - line ' ...
            num2str(stack(1).line) ...
            ': Starting graph analysis ...'];
    ft_progress('init', 'text', msg)
end

% count number of cases were either method doesn't return an alternative path
n_nopath_TEdyn          = 0;
n_nopath_TEbacktracking = 0;

for i=1:n_edges;
    if ~strcmp(cfg.verbosity, 'none')
        ft_progress(i/n_edges, [repmat('   ', 1, length(dbstack)-1) '   Processing edge %d of %d ...'], i, n_edges);
    end
    
    % define current source, target and upper limit k
    k = weights(i) + threshold;      
    s = edges4analysis(i,1);
    t = edges4analysis(i,2);    
    
    if k<=0;
        error('Something is wrong with your threshold!');
    end;
    
    % remove current edge for this run of the algorithm
    edges_temp = edges4analysis;
    edges_temp(i,:) = [];
    weights_temp = weights;
    weights_temp(i) = [];
    
    % rearrange labels, thus source=1 and target=end, the mapping, enumeration can be changed back after
    % backtracking  
    labels_vertices_temp = 1:n_vertices;
    labels_vertices_temp(labels_vertices_temp==s) = [];
    labels_vertices_temp(labels_vertices_temp==t) = [];
    labels_vertices_temp = [s; labels_vertices_temp'; t];
    % mask is needed for the new enumeration of the vertices
    mask    = ones(size(edges_temp));
    % 
    
    % enumerate nodes, masking is needed so that already changed nodes,
    % don't get changed again (happens if a node is changed to a higher
    % number n and if j=n, this node is overwritten again)
    for j=1:length(labels_vertices_temp);
        mask_temp = edges_temp==labels_vertices_temp(j)&mask;
        edges_temp(mask_temp) = j;
        mask(mask_temp) = 0;
    end;
    clear mask mask_temp;
        
    % create 'inverted' adjacency list
    adjacency_list = cell(n_vertices,1);
    for j=1:n_vertices;
        % check whether node has any predecessors
        if ~isempty(edges_temp(edges_temp(:,2)==j,1))
            % add predecessors for current node j
            adjacency_list{j} = cat(1, ...
                edges_temp(edges_temp(:,2)==j,1)', ...      % find all predecessors of vertex j
                weights_temp(edges_temp(:,2)==j)');         % find corresponding edge-weights
        end;
    end;

	% if s and t are part of the same subgraph, look for alternative paths
	if(TEdfs(adjacency_list))
        
		solution = TEdyn(adjacency_list,k);
        
        % check if alternative paths were found
        alt_paths = 0;
        for j=k-2*threshold:k;
            if j<1; continue; end;
            if ~isempty(solution{j+1,end});
                alt_paths = 1;
                break;
            end;
        end;
        
        % if alternative paths exist, do backtracking
        if logical(alt_paths);
            
            paths = TEbacktracking(solution,k,threshold);

            if ~isempty(paths)   

                % change enumeration back to original format (before
                % deletion of current edge) and count alternative paths
                path_count = size(paths,2);               
                for j=1:path_count; 
                    paths{j} = labels_vertices_temp(paths{j});                    
                end;
                
                % collect alternative paths
                all_paths = [all_paths; {i s t path_count paths}];   

            else
                %disp('  No alternative paths found by TEbacktracking.');
                n_nopath_TEbacktracking = n_nopath_TEbacktracking+1;
            end;
        else
            %disp('  No alternative paths found by TEdyn.');
            n_nopath_TEdyn          = n_nopath_TEdyn+1;            
        end;
    %else
        %disp('  Source and target are not in the same subgraph');
 	end;

end;
if ~strcmp(cfg.verbosity, 'none')
    ft_progress('close');
end
%% prepare output

% flag all edges to which alternative paths exist
if ~isempty(all_paths)
    
    msg = sprintf('Alternative paths were found for %d of %d edges', size(all_paths,1), n_edges);
    TEconsoleoutput(cfg.verbosity, msg, LOG_INFO_MINOR);
    
    % add alternative paths and graph info to datastructure
    [data_paths, triangle_edges, triangle_nodes] = TEflagedges(data,all_paths,edges4analysis,edges_original);    
    data_paths.graphanalysis = graphanalysis;
    data_paths.graphanalysis.cfg = cfg;
    data_paths.graphanalysis.triangle_edges = triangle_edges;
    data_paths.graphanalysis.triangle_nodes = triangle_nodes;
    
    msg = sprintf('%d triangle(s) were found by TEflagedges', size(triangle_edges, 1));
    TEconsoleoutput(cfg.verbosity, msg, LOG_INFO_MINOR);
    
else
    data_paths = data;
    data_paths.graphanalysis = graphanalysis;
    TEconsoleoutput(cfg.verbosity, 'No alternative paths were found!', LOG_INFO_MINOR);
end;

% update TEsteps
if ~isfield(data,'TEsteps')      %adding structure with changings; added modified by nicu
    data_paths.TEsteps = 'GA';
else
    data_paths.TEsteps = strcat(data.TEsteps,'_GA');
end

ft_progress('close');
