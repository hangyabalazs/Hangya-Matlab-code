function [Burstiness,IntraBurstFrequency,IntraBurstSpikeNumber,BurstLength,BurstFrequency,Silence,...
    Gap,After,IsiMatrix,Burst,H] = aclust(vdisc,method)
%ACLUST   Hierarchical clustering.
%   ACLUST(VDISC,METHOD) hierarchically clusters the input argument VDISC (discriminated unit). 
%   The clustering method is set to 'ward' by default and can be set to 'average' or 'single' 
%   using METHOD parameter. See header line of the program code for detalied call of the function.
%
%   Outputs:
%       1. Burst parameters (burstiness (number of intraburst spikes / number of all spikes), intraburst
%           frequency (calculated for each burst and averaged) - mean, sd and separately for all bursts,
%           intraburst spike number - mean, sd and separately for all bursts, burst length - mean, sd 
%           and separately for all bursts, burst frequency).
%       2. Silence (minimum preburst silent periods), gap (silence minus maximal intraburst interval),
%           after (afterburst periods)
%       3. IsiMatrix (first, second, etc. interspike intervals in bursts with different spike number
%           - mean, sd and number of bursts with given spike number)
%       4. Burst cell (a cell array containing the burst limits for every structural element (substate).
%       5. H (figure handle; plot contains intraburst frequency, silence, gap and reciprocial of
%           maximal burst-first-iv. vs number of clusters)
%
%   See also ITCLUST2, CLUSTER, COPHENET, INCONSISTENT and NATCLUST.

% Input arguments check
error(nargchk(1,2,nargin));
if nargin == 1
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

% Allocating output variables
Burst = cell(1,livs);
After = cell(1,livs);
Silence = zeros(1,livs);
Gap = zeros(1,livs);
Burstiness = zeros(1,livs);
IntraBurstFrequency = struct('mean',[],'sd',[],'all',{});
ibf_mean = zeros(1,livs);
ibf_min = zeros(1,livs);
firstiv_max = zeros(1,livs);
IntraBurstSpikeNumber = struct('mean',[],'sd',[],'all',{});
BurstLength = struct('mean',[],'sd',[],'all',{});
BurstFrequency = zeros(1,livs);
IsiMatrix = struct('mean',[],'sd',[],'num',[]);

% CLUSTER
wb = waitbar(0,'Running ACLUST for Andi to be happy :-)...','Position',[360 250 275 50]);  %progress indicator
global WB
WB(end+1) = wb;

miniv = min(ivs);   %finding the cluster containing the smallest interval
for dec = 2:livs
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
    firstiv = [];
    burstnum = size(burst,2);
    intraburstnum = zeros(1,burstnum);
    isimtx = {};
    for j = 1:burstnum    %computing intraburstiv and allfirstspike
        b = vdisc(burst(1,j):burst(2,j));
        db = diff(b);
        firstiv = [firstiv db(1)];
        intraburstiv = [intraburstiv db];                 
        ssi(burst(1,j)+1:burst(2,j)) = 0;
        intraburstnum(j) = length(b);   %intraburst spike number
        for k = 1:length(b) - 1     % ISI matrix
            if size(isimtx,1) < (length(b) - 1) | size(isimtx,2) < k
                isimtx{length(b)-1,k} = [];
            end
            isimtx{length(b)-1,k} = [isimtx{length(b)-1,k} (vdisc(burst(1,j)+k)-vdisc(burst(1,j)+k-1))/20000];
        end
    end
    burstlength = (vdisc(burst(2,:)) - vdisc(burst(1,:))) / 20000;
    intraburstfreq = (intraburstnum - 1) ./ burstlength;
    
% Interburst intervals
    fsp2 = fsp(2:end);  %computing interburstiv
    lsp2 = lsp(1:end-1);
    interburstiv = fsp2 - lsp2;
    
% Inter-1st-sike intervals
    interfirstspike = diff(fsp);
    if ~isempty(interfirstspike)
        firstspikefreq = 20000 * (burstnum - 1)  / (vdisc(burst(2,end)) - vdisc(burst(1,1)));
    else
        firstspikefreq = NaN;
    end
    
% Beforefirst intervals
    if ~isequal(burst(1,1),1)
        beforefirstiv = ivs(burst(1,:)-1);
    else
        beforefirstiv = ivs(burst(1,2:end)-1);
    end
    Silence(dec) = min(beforefirstiv) / 20000;
    gap = (min(beforefirstiv) - max(intraburstiv)) / 20000;
    Gap(dec) = gap;
    
% Afterlast intervals
    if ~isequal(burst(end,end),length(vdisc))
        afterlastiv = ivs(burst(2,:));
    else
        afterlastiv = ivs(burst(2,1:end-1));
    end
    After{dec} = afterlastiv / 20000;
    
% Burst parameters
    Burstiness(dec) = (length(intraburstiv) + burstnum) / length(vdisc);
    IntraBurstFrequency(dec).mean = mean(intraburstfreq);
    ibf_mean(dec) = mean(intraburstfreq);
    ibf_min(dec) = min(intraburstfreq);
    firstiv_max(dec) = max(firstiv);
    IntraBurstFrequency(dec).sd = std(intraburstfreq);
    IntraBurstFrequency(dec).all = intraburstfreq;
    IntraBurstSpikeNumber(dec).mean = mean(intraburstnum);
    IntraBurstSpikeNumber(dec).sd = std(intraburstnum);
    IntraBurstSpikeNumber(dec).all = intraburstnum;
    BurstLength(dec).mean = mean(burstlength);
    BurstLength(dec).sd = std(burstlength);
    BurstLength(dec).all = burstlength;
    BurstFrequency(dec) = firstspikefreq;
    for x = 1:size(isimtx,1)
        warning off
        for y = 1:size(isimtx,2)
            IsiMatrix(dec).mean(x,y) = mean(isimtx{x,y});
            IsiMatrix(dec).sd(x,y) = std(isimtx{x,y});
            IsiMatrix(dec).num(x,y) = length(isimtx{x,y});
        end
        warning backtrace
    end
    
    waitbar(dec/livs)
end
close(wb)

% Plot
H = figure;
subplot(4,1,1)
plot(ibf_mean)
subplot(4,1,2)
plot(Silence)
subplot(4,1,3)
plot(Gap)
subplot(4,1,4)
plot(ibf_min)

% If no outputs are expected uninitialize all outputs to prevent ans displays
if nargout < 1
    clear
end