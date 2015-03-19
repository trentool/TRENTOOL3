function solution = TEdyn(adjacency_list,k);
% TEDYN: Dynamic programming algorithm to find all paths up to a given 
% length k. The first node in the adjaency list is assumed to be the 
% source s, the last node to be the target t of all paths.
%
% * INPUT PARAMETERS
%       adjacency_list = predecessor structure: 'inverted' adjacency list,
%                        where for every node all predecessors are listed
%       k              = upper bound of reconstruction interval
%                        (w_crit+threshold)
%
%
% * OUTPUT PARAMETERS
%       solution       = array that contains recurrent solutions to all
%                        subproblems of the problem whether a path of
%                        length w_crit +/- a threshold exists
%
% See the reference information for a more detailed description.
%
%
% * REFERENCE INFORMATION
%
%   - graph algorithm
%         - Bsc Thesis Patricia Wollstadt 
%           (email: p.wollstadt@stud.uni-frankfurt.de)
%   - dynamic programming
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
% PW - 06/09/2011

n_vertices = size(adjacency_list,1);

solution = cell(k+1,n_vertices); % initialize solution-array
solution{1,1}(1) = 0; % s is reachable by a path of length 0

% loop over rows and columns of solution-array
for i=1:k;
    for v=1:n_vertices;
        %disp(['Solving subproblem k=' num2str(i) '/v=' num2str(v)]);
        % loop over all edges pointing to current vertex
        for e=1:size(adjacency_list{v},2);
            u = adjacency_list{v}(1,e); % define current parent
            w = round(adjacency_list{v}(2,e));	% define corresponding edge-weight
            if (i-w) >= 0 
                try
                if (~isempty(solution{i-w+1,u}))
                    solution{i+1,v} = cat(2,solution{i+1,v},[u;i-w]); % remember current parent, if path to parent exists
                end;
                catch
                    disp('error')
                end; 
            end;
        end;
    end;
end;
