function [paths_return] = TEbacktracking(solution_TEdyn,k,threshold)
% TEBACKTRAKING:Backtracking for solutions produced by TEdyn. Function 
% reconstructs all paths to from source to target node with lengths within 
% the given threshold (w_crit +/- the user-defined threshold theta).
%
%
% * INPUT PARAMETERS 
%       solution_TEdyn = solution array returned by TEdyn
%       k              = upper bound of reconstruction interval
%                        (w_crit+threshold)
%       threshold      = threshold around w_crit, allows for imprecision in
%                        the reconstruction of delay times
%
% * OUTPUT PARAMETERS
%       paths           = cell array, that contains cell arrays of
%                         alternative paths 
%
%
% Refer to the reference information for a more detailed description.
%
% * DEPENDENCIES
%       - functions
%           - TEbacktracking_rec (part of this file, see below)
%
%
% * REFERENCE INFORMATION
%
%   - graph algorithm
%         - Bsc Thesis Patricia Wollstadt 
%           (email: p.wollstadt@stud.uni-frankfurt.de)
%
% PW - 06/09/2011
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation;
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY;
%
% Version 1.2 by Patricia Wollstadt
% Frankfurt 2012

% CHANGELOG
% 06-23-14: PW changed the output format


solution_tree = cell(2*threshold+1,2);
count=0;

% check if lower limit is below zero (happens with small edge-weights) 
lower_limit = k-2*threshold;
if lower_limit<1;
    lower_limit = 1;
    %disp('  Lower limit was set to 1.');
else
    %disp(['  Lower limit is ' num2str(lower_limit)]);
end;

% initialize visited list
visited = zeros(1,size(solution_TEdyn,2));
visited(1) = 1; % set source to visited

paths_return = {};

% loop over all path-lengths within threshold
for i=lower_limit:k;	
	s = [i+1 size(solution_TEdyn,2)]; % determine current source, i+1 accounts for matlab-indexing 
    if ~isempty(solution_TEdyn{s(1),s(2)}) % if solution exists, do backtracking        
        paths = TEbacktracking_rec(s,solution_TEdyn,{},[],visited,0);
        if ~isempty(paths);

            for j = 1:length(paths); % sort paths and add source
                paths{j} = cat(2,paths{j},1);
                paths_return = cat(2,paths_return,fliplr(paths{j}));
            end;
            count = count+1;
            solution_tree{count,1} = paths; % add paths to result
            solution_tree{count,2} = i;     % add path-weight to result
        end;
    end;
end;

% just keep non-empty cells
%paths     = solution_tree(1:count,1);
%paths_len = [solution_tree{1:count,2}];            % PW: this is currently not needed in any subsequent analysis step - comment in if needed

% if isempty(solution_tree)
%     %disp('No Paths found.');
%     paths     = [];
%     %paths_len = [];
% end;



function [paths current_path] = TEbacktracking_rec(source,solution,paths,current_path,visited,depth)
% function [paths current_path] = TEbacktracking_rec(source,solution,paths,current_path,visited)
%
% Uses a recursive depth-first approach to backtrack all paths from target 
% to source. Function is exclusively called by 
% TEbacktracking(solution_TEdyn,k,threshold).
%
% PW - 06/09/2011

% base cases
if (source==[1 1])  % case 1: starting node is reached by a path of length k
    %disp('Valid path found.');
	paths = cat(1,paths,current_path);
	return;
elseif (visited(source(2)) == 1) % case 2: current node has been visited before, path contains a loop 
    %disp('Path contains a loop. Return.');
    return;
elseif length(paths) > 20000 % case 3: too many paths
    disp('Too many paths. Return.');
    return;
end;

% recursive call, if no base case
for i=1:size(solution{source(1),source(2)},2)
    
    % add current node to current path and mark as visited
    current_path = [current_path source(2)];
    visited(source(2)) = 1;
    %disp(['Recursive depth: ' num2str(depth) ' Current Path ' num2str(current_path)]);
    depth = depth + 1; 
    new_source = [solution{source(1),source(2)}(2,i)+1 solution{source(1),source(2)}(1,i)];
	[paths current_path] = TEbacktracking_rec(new_source,solution,paths,current_path,visited,depth);
    depth = depth-1;
    visited(source) = 0;
    current_path(end) = []; 
end;