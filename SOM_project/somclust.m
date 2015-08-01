function [c cc] = somclust(fcns,dec,gr1,gr2,method)
%SOMCLUST   Hierarchical clustering.
%   C = SOMCLUST(V,DEC,G1,G2,METHOD) hierarchically clusters the input
%   argument V into n = DEC clusters (default = 5). The clustering method
%   is set to 'ward' by default and can be set to 'average' or 'single'
%   using METHOD parameter. G1 and G2 define groups for which leaves and
%   branches appear in color (see DENDROGRAM_COLORED).
%
%   See also ACLUST, DENDROGRAM_COLORED, ITCLUST2, CLUSTER, COPHENET,
%   INCONSISTENT and NATCLUST.

% Input arguments check
error(nargchk(1,5,nargin));
if nargin < 5
    method = 'ward';     %set default method to 'ward';
end
if nargin < 2
    dec = 5;    % default number of clusters
end

% DIST, LINKAGE and COPHENET
livs = length(fcns);
dist = pdist2(fcns);
links = linkage(dist,method);
cc = cophenet(links,dist);

% CLUSTER
c = cluster(links,dec);
% links(:,3) = log(links(:,3));   % log of distances for the dendrogram
if nargin > 2
    dendrogram_colored(links,gr1,gr2,0);
else
    dendrogram(links,0);
end