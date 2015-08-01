function [intraburstiv_cv,extraburstiv_cv,interburstiv_cv,...
        firstspike_cv,allfirstspike_cv,burstlength_cv,...
        Burstiness,IntraBurstFrequency,IntraBurstSpikeNumber,...
        BurstLength,BurstFrequency,Burst,Ccc] = b_itclust2(vdisc,d,method)
%ITCLUST2   Iterative clustering.
%   ITCLUST2(VDISC,D,METHOD) clusters the input argument VDISC (discriminated unit) using the "iterative 
%   clustering strategy" setting the number of clusters to D. The clustering method is set to 'ward' 
%   by default and can be set to 'average' or 'single' using METHOD parameter. See first line of the
%   program code for detalied call of the function.
%
%   Outputs:
%       1. Coefficient of variation (CV) values (intraburst interval, iterburst interval, extraburst
%           interval, "inter-first-spike interval" ,"all-inter-first-spike-interval" (distances of
%           first spikes including single spikes as well), burst length CV).
%       2. Burst parameters (burstiness (number of intraburst spikes / number of all spikes), average
%           intraburst frequency (calculated for each burst and averaged), average intraburst spike.
%           number, average burst length, frequency of burst first spikes (called 'burst frequency').
%       3. Burst cell (a cell array containing the burst limits for every structural element (substate).
%       4. Cophenetic coefficient ('ccc' - see COPHENET for details).
%
%   See also ICA_NEWAVE, CLUSTER, COPHENET, INCONSISTENT and NATCLUST.

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

% Allocating some matrices
intraburstiv_cv = zeros(1,d);  
extraburstiv_cv = zeros(1,d);
interburstiv_cv = zeros(1,d);
firstspike_cv = zeros(1,d);
allfirstspike_cv = zeros(1,d);

Burst = cell(1,d);  %allocating Burst cell

% CLUSTER
miniv = min(ivs);   %finding the cluster containing the smallest interval
for dec = 2:d
    c = cluster(links,dec);
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
    
    Burst{dec} = burst; %filling Burst cell
    
% Intraburst intervals
    ssi = vdisc;    %ssi is for "single spikes included": it will contain the first spikes of bursts
                    %and the single spikes as well
    intraburstiv = [];
    sb2 = size(burst,2);
    intraburstnum = zeros(1,sb2);
    for j = 1:sb2    %computing intraburstiv and allfirstspike
        b = vdisc(burst(1,j):burst(2,j));
        intraburstiv = [intraburstiv diff(b)];                 
        ssi(burst(1,j)+1:burst(2,j)) = 0;
        intraburstnum(j) = length(b);   %intraburst spike number
    end
    burstlength = (vdisc(burst(2,:)) - vdisc(burst(1,:))) / 10000;
    intraburstfreq = (intraburstnum - 1) ./ burstlength;
    
% Interburst intervals
    fsp2 = fsp(2:end);  %computing interburstiv
    lsp2 = lsp(1:end-1);
    interburstiv = fsp2 - lsp2;
    
% Inter-1st-sike intervals
    interfirstspike = diff(fsp);
    if ~isempty(interfirstspike)
        firstspikefreq = 10000 / mean(interfirstspike);
    else
        firstspikefreq = NaN;
    end
    
% CV of the intraburst interval length
    if length(intraburstiv) ~= 0
        intraburstiv_cv(dec) = std(intraburstiv) / mean(intraburstiv);  
    else intraburstiv_cv(dec) = NaN;
    end
    %an alternative way of calculating intraburstiv_cv:
    %intraburstiv_cv2 = std(ivss)/mean(ivss);
    
% CV of the extraburst interval length
    if length(extraburstiv) ~= 0
        extraburstiv_cv(dec) = std(extraburstiv) / mean(extraburstiv);
    else extraburstiv_cv(dec) = NaN;
    end
    
% CV of the interburst interval length
    if length(interburstiv) ~= 0
        interburstiv_cv(dec) = std(interburstiv) / mean(interburstiv);
    else interburstiv_cv(dec) = NaN;
    end
    
% CV of the inter-1st-spike interval length
    if length(interfirstspike) ~= 0
        firstspike_cv(dec) = std(interfirstspike) / mean(interfirstspike);
    else firstspike_cv(dec) = NaN;
    end
    
% CV of the all-inter-1st-spike interval length
    dfssi = diff(ssi(find(ssi)));
    allfirstspike_cv(dec) = std(dfssi) / mean(dfssi);    %Note, that dfssi cannot be empty!
    
% CV of burst length
    if length(burstlength) ~= 0
        burstlength_cv(dec) = std(burstlength) / mean(burstlength);
    else burstlength_cv(dec) = NaN;
    end

% Burst parameters
    Burstiness(dec) = (length(intraburstiv) + size(burst,2)) / length(vdisc);
    IntraBurstFrequency(dec) = mean(intraburstfreq);
    IntraBurstSpikeNumber(dec) = mean(intraburstnum);
    BurstLength(dec) = mean(burstlength);
    BurstFrequency(dec) = firstspikefreq;
end

% If no outputs are expected uninitialize all outputs to prevent ans displays
if nargout < 1
    clear
end