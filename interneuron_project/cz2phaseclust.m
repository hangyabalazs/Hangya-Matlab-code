function Ccc = cz2phaseclust(angs,method)
%CZ2PHASECLUST   Clustering for phase values.
%   CCC = CZ2PHASECLUST(ANGS,METHOD) clusters the circular input argument ANGS.
%   It uses the "natural clustering strategy" setting the cutoff inconsistency coefficient to different ICC values. The clustering method
%   is set to 'average' by default. 
%
%   Output: Cophenetic coefficient ('ccc' - see COPHENET for details).
%
%   See also CLUSTER, COPHENET, INCONSISTENT and NATCLUST.

% Input arguments check
error(nargchk(2,3,nargin));
if nargin == 2
    method = 'average';    %set default method to 'ward';
end

% DIST, LINKAGE and COPHENET
ivs = angs;
livs = length(ivs);
dmtx = zeros(livs,2);
iivs = ivs';
dmtx(:,1) = iivs;
dist = pdist2(dmtx,'circular');
links = linkage(dist,method);
Ccc = cophenet(links,dist);

c = cluster(links,'maxclust',4);
figure
plot(c,angs*180/pi,'.')

% CLUSTER
clusno = [];
for ICC = 0.5:0.05:2
    c = cluster(links,'cutoff',ICC);
    dec = max(c);
    clusno(end+1) = dec;    %number of clusters
end
figure
plot(clusno)