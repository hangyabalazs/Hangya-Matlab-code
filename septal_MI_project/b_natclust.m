function [intraburstiv_cv,extraburstiv_cv,interburstiv_cv,firstspike_cv,allfirstspike_cv,...
        Burstiness,IntraBurstFrequency,IntraBurstSpikeNumber,BurstLength,BurstFrequency,...
        Burst,Ccc,clusno] = b_natclust(vdisc,ICC,method)
%NATCLUST   Natural clustering.
%   [BURST,CCC] = NATCLUST(VDISC,ICC,METHOD) clusters the input argument VDISC using the "natural 
%   clustering strategy" setting the cutoff inconsistency coefficient to ICC. The clustering method
%   is set to 'ward' by default and can be set to 'average' or 'single' using METHOD parameter. 
%
%   Outputs:
%       1. Coefficient of variation (CV) values (intraburst interval, iterburst interval, extraburst
%           interval, "inter-first-spike interval" and "all-inter-first-spike-interval" (distances of
%           first spikes including single spikes as well) CV).
%       2. Burst parameters (burstiness (number of intraburst spikes / number of all spikes), average
%           intraburst frequency (calculated for each burst and averaged), average intraburst spike.
%           number, average burst length, frequency of burst first spikes (called 'burst frequency').
%       3. Burst cell (a cell array containing the burst limits for every structural element (substate).
%       4. Cophenetic coefficient ('ccc' - see COPHENET for details).
%       5. Number of created clusters.
%
%   See also ICA_NEWAVE, CLUSTER, COPHENET and INCONSISTENT.

% Input arguments check
error(nargchk(2,3,nargin));
if nargin == 2
    method = 'ward';    %set default method to 'ward';
end

% DIST, LINKAGE and COPHENET
ivs = diff(vdisc);
livs = length(ivs);
dmtx = zeros(livs,2);
iivs = ivs';
dmtx(:,1) = iivs;
dist = pdist2(dmtx);
links = linkage(dist,method);
Ccc = cophenet(links,dist);

% CLUSTER
miniv = min(ivs);   %finding the cluster containing the smallest interval
c = cluster(links,ICC);
dec = max(c);
clusno = dec;    %number of clusters

% Finding the cluster containing the smallest interval
cell_clusters = cell(dec,1);
doub_length = zeros(dec,1);
for t = 1:dec
    cell_clusters{t} = find(c==t);
    doub_length(t) = length(cell_clusters{t});
    fnd = find(ivs(cell_clusters{t}) == miniv);
    if isempty(fnd) == 0,
        miniv_clus = t;
    end
end

% Extraburst intervals
ivss = ivs(cell_clusters{miniv_clus});  %creating extraburstiv
extraburstiv = ivs(find(ivs>max(ivss)));

% Bursts
liv = find(ivs<=max(ivss)); %finding the bursts
liv1 = [-1 liv]; 
liv = [liv 0];
bliv = find(liv~=liv1+1);
blivv = bliv(2:end)-1; 
bliv(end) = [];
burst = [liv(bliv);liv(blivv)+1];     %1st and last spikes of bursts
fsp = vdisc(burst(1,:));              %localisation of 1st spikes of bursts
lsp = vdisc(burst(2,:));              %localisation of last spikes of bursts
diffburst = vdisc(burst(2,:)) - vdisc(burst(1,:));
mdb = max(diffburst);

Burst = burst;

% Intraburst intervals
ssi = vdisc;    %ssi is for "single spikes included": it will contain the first spikes of bursts
                %and the single spikes as well
intraburstiv = [];
for j = 1:size(burst,2),    %computing intraburstiv and allfirstspike
    b = vdisc(burst(1,j):burst(2,j));
    intraburstiv = [intraburstiv diff(b)];                   
    ssi(burst(1,j)+1:burst(2,j)) = 0;
    intraburstnum(j) = length(b);   %intraburst spike number
end
burstlength = vdisc(burst(2,:)) - vdisc(burst(1,:));
intraburstfreq = (intraburstnum - 1) ./ burstlength;

% Interburst intervals
fsp2 = fsp(2:end);  %computing interburstiv
lsp2 = lsp(1:end-1);
interburstiv = fsp2 - lsp2;

% Inter-1st-sike intervals
interfirstspike = diff(fsp);
firstspikefreq = 1 / mean(interfirstspike);

% CV of the intraburst interval length
if length(intraburstiv) ~= 0
    intraburstiv_cv = std(intraburstiv)/mean(intraburstiv);  
else intraburstiv_cv = NaN;
end
%an alternative way of calculating intraburstiv_cv:
%intraburstiv_cv2 = std(ivss)/mean(ivss);

% CV of the extraburst interval length
if length(extraburstiv) ~= 0
    extraburstiv_cv = std(extraburstiv)/mean(extraburstiv);
else extraburstiv_cv = NaN;
end

% CV of the interburst interval length
if length(interburstiv) ~= 0
    interburstiv_cv = std(interburstiv)/mean(interburstiv);
else interburstiv_cv = NaN;
end

% CV of the inter-1st-spike interval length
if length(interfirstspike) ~= 0
    firstspike_cv(dec) = std(interfirstspike) / mean(interfirstspike);
else firstspike_cv(dec) = NaN;
end

% CV of the all-inter-1st-spike interval length
dfssi = diff(ssi(find(ssi)));
allfirstspike_cv = std(dfssi)/mean(dfssi);    %Note, that dfssi cannot be empty!

% Burst parameters
Burstiness(dec) = (length(intraburstiv) + size(burst,2)) / length(vdisc);
IntraBurstFrequency(dec) = mean(intraburstfreq);
IntraBurstSpikeNumber(dec) = mean(intraburstnum);
BurstLength(dec) = mean(burstlength);
BurstFrequency(dec) = firstspikefreq;

% If no outputs are expected uninitialize all outputs to prevent ans displays
if nargout < 1
     clear
end