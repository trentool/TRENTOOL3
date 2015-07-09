function dens = getDensity(V,E)
% GETDENSITY: Calculates the density of a graph as the ratio of the
% existing and potential edges in a directed graph (weighting is not taken 
% into account).
%
% * INPUT
%    EITHER
%       V = no. vertices AND 
%       E = no. edges
%   OR
%       adjacency matrix = matrix that represents the TE graph, where the
%                          field in the i-th row and j-th column is the
%                          weight of the edge from the i-th to the j-th 
%                          vertex
%
% * OUTPUT
%       dens = density of the graph as described above
%
% * REFERENCE INFORMATION
%       - S. Wasserman (1994). Social Network Analysis: Methods and 
%         Applications. Cambridge University Press.
%
% PW 16/10/11
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

if nargin == 1;   
    E = sum(sum(V>0)); 
    V = size(V,1);
end;

dens = E/(V*(V-1));