function bool_connected = TEdfs(adjacency_list)
% function bool_connected = TEdfs(adjacency_list)
%
% Checks if source (last vertex in adjacency list) and target (first
% vertex in adjacency list) are part of the same subgraph. For this a depth
% first search (DFS) is used.
%
% * INPUT PARAMETERS
%       adjacency_list = cell array (no. vertices x 1) of vectors, where 
%                        each entry in the cell array is treated as a node
%                        and the vector is a list of the children of this
%                        node (represents the complete TE graph)
%
% * OUTPUT PARAMETERS
%       bool_connected = boolean indicating, whether the first and last
%                        node are connected (1) or not (0)
%
%
% * REFERENCE INFORMATION
%
%   - depth first search
%         - the basic principle is for example described in Cormen et. al. 
%           (2009). Introduction to Algorithms (third edition). The MIT 
%           Press.
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
%
% PW - 07/09/2011

source = length(adjacency_list); 
target = 1;

bool = TEdfs_rec(source,target,adjacency_list,zeros(1,length(adjacency_list)),0);

bool_connected = logical(bool);



function bool = TEdfs_rec(source,target,adjacency_list,visited,bool)
% function bool = TEdfs_rec(source,target,adjacency_list,visited,bool)
%
% Uses a recursive depth-first approach to check, if s and t are in the
% same subgraph. Function is exclusively called by TEdfs(adjacency_list).
%
% PW - 07/09/2011

if source == target; % base case: target is reached
    bool = 1;
    return;
else
    for i=1:size(adjacency_list{source},2); % check all children of current source
        if visited(adjacency_list{source}(1,i)) == 0; % recursive call, if current child has not been visited yet
            new_source = adjacency_list{source}(1,i);
            visited(new_source) = 1;            
            bool = TEdfs_rec(new_source,target,adjacency_list,visited,bool);
        end;
        if logical(bool); return; end; % abort recursion, if target is reached
    end;
end;

