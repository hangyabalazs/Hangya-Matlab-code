function somentropy_optimize_wn
%SOMENTROY_OPTIMIZE_WN   Optimize cutoff values for spike width.
%   SOMENTROPY_OPTIMIZE_WN optimizes cutoff values of spike width for
%   narrow spiking and wide spiking cells (separately) for information
%   calculations (see SOMENTROPY_OPTIMIZE2B) over the range of 180-380 us.
%
%   See also SOMENTROPY_OPTIMIZE2B.

% Load & define PSTHs
global DATADIR
load([DATADIR 'SOM_Sachin\PSTH\HomeZout.mat'])
try   % different naming conventions
    pv_MatrixPsth1 = pv_psth;
    som_MatrixPsth1 = som_psth;
    non_tagged_MatrixPsth1 = non_tagged_all;
catch
    pv_MatrixPsth1 = pv_psth;
    som_MatrixPsth1 = som_psth;
    non_tagged_MatrixPsth1 = nt_psth;
end

load([DATADIR 'SOM_Sachin\PSTH\HomeZin.mat'])
try   % different naming conventions
    pv_MatrixPsth2 = pv_psth;
    som_MatrixPsth2 = som_psth;
    non_tagged_MatrixPsth2 = non_tagged_all;
catch
    pv_MatrixPsth2 = pv_psth;
    som_MatrixPsth2 = som_psth;
    non_tagged_MatrixPsth2 = nt_psth;
end

% Load spike width data
load([DATADIR 'SOM_Sachin\PSTH\SpikeWidth.mat'])
if ~exist('PvSpikeWidth','var')
    NonTaggedAllSpikeWidth = nt_cells;
    PvSpikeWidth = pv;
    SomSpikeWidth = som;
end

% Remove NaNs
[x1 y1] = find(isnan(non_tagged_MatrixPsth1));
[x2 y2] = find(isnan(non_tagged_MatrixPsth2));
x = union(unique(x1),unique(x2));
non_tagged_MatrixPsth1(x,:) = [];
non_tagged_MatrixPsth2(x,:) = [];
NonTaggedAllSpikeWidth(x) = [];
x = find(isnan(NonTaggedAllSpikeWidth));
non_tagged_MatrixPsth1(x,:) = [];
non_tagged_MatrixPsth2(x,:) = [];
NonTaggedAllSpikeWidth(x) = [];

% Restrict in time
wn1 = -60;
wn2 = 60;
mdl = (size(pv_MatrixPsth1,2) + 1) / 2;
pv_MatrixPsth1 = pv_MatrixPsth1(:,mdl+wn1:mdl+wn2-1);
som_MatrixPsth1 = som_MatrixPsth1(:,mdl+wn1:mdl+wn2-1);
non_tagged_MatrixPsth1 = non_tagged_MatrixPsth1(:,mdl+wn1:mdl+wn2-1);
pv_MatrixPsth2 = pv_MatrixPsth2(:,mdl+wn1:mdl+wn2-1);
som_MatrixPsth2 = som_MatrixPsth2(:,mdl+wn1:mdl+wn2-1);
non_tagged_MatrixPsth2 = non_tagged_MatrixPsth2(:,mdl+wn1:mdl+wn2-1);

% Smooth
sizepv = size(pv_MatrixPsth1);
sizesom = size(som_MatrixPsth1);
sizent = size(non_tagged_MatrixPsth1);
spv1 = smoothallcells(pv_MatrixPsth1);   % PV, HomeZoneOut
ssom1 = smoothallcells(som_MatrixPsth1);   % SOM, HomeZoneOut
snt1 = smoothallcells(non_tagged_MatrixPsth1);   % non-tagged, HomeZoneOut
spv2 = smoothallcells(pv_MatrixPsth2);   % PV, HomeZoneIn
ssom2 = smoothallcells(som_MatrixPsth2);   % SOM, HomeZoneIn
snt2 = smoothallcells(non_tagged_MatrixPsth2);   % non-tagged, HomeZoneIn

% HomeZoneOut
som_merged = ssom1;
pv_merged = spv1;
nt_merged = snt1;
allpsth = [pv_merged; som_merged; nt_merged];

% Convert PSTHs to numbers
[allpsth5_pv nbins_pv] = convertpsth(allpsth,pv_merged);

% Tag distribution
tags_pv = [true(size(pv_merged,1),1); false(size(som_merged,1)+size(nt_merged,1),1)];  % knowing PV

% Kullback-Leibler divergence
klb_pv = KLdist_discrete(allpsth5_pv(tags_pv),allpsth5_pv);

% Correction with surrogate
sss = 100;    % size of surrogate sample
msrg = calc_surrogate(allpsth,sss,tags_pv,nbins_pv);
kl_pv = klb_pv - msrg;

% Narrow spiking cells
L = 180:10:360;
kl_ns = nan(size(L));
klb_ns = nan(size(L));
next = 1;
for lmt = L
    disp(lmt)
    
    % Tag distribution
    tags_ns = [PvSpikeWidth<lmt; SomSpikeWidth<lmt; NonTaggedAllSpikeWidth<lmt];  % knowing narrow spiking
    
    % Correlate with mean PV PSTH
    [allpsth5_ns nbins_ns] = convertpsth(allpsth,allpsth(tags_ns,:),nbins_pv);   % note: nbins_pv is used for NS and WS as well
    
    % Kullback-Leibler divergence
    klb_ns(next) = KLdist_discrete(allpsth5_ns(tags_ns),allpsth5_ns);
    
    % Correction with surrogate
    msrg = calc_surrogate(allpsth,sss,tags_ns,nbins_ns);
    kl_ns(next) = klb_ns(next) - msrg;
    
    next = next + 1;
end

% Wide spiking cells
kl_ws = nan(size(L));
klb_ws = nan(size(L));
next = 1;
for lmt = L
    disp(lmt)
    
    % Tag distribution
    tags_ws = [PvSpikeWidth>lmt; SomSpikeWidth>lmt; NonTaggedAllSpikeWidth>lmt];  % knowing wide spiking
    
    % Correlate with mean PV PSTH
    [allpsth5_ws nbins_ws] = convertpsth(allpsth,allpsth(tags_ws,:),nbins_pv);
    
    % Kullback-Leibler divergence
    klb_ws(next) = KLdist_discrete(allpsth5_ws(tags_ws),allpsth5_ws);
    
    % Correction with surrogate
    msrg = calc_surrogate(allpsth,sss,tags_ws,nbins_ws);
    kl_ws(next) = klb_ws(next) - msrg;
    
    next = next + 1;
end

% Optimize cut for spike width
sumkl = kl_ns + kl_ws;
figure;
plot(L,sumkl)
inx = sumkl == max(sumkl);
lmt = L(inx);
kl_nso = kl_ns(inx);
kl_wso = kl_ws(inx);
klb_nso = klb_ns(inx);
klb_wso = klb_ws(inx);
tags_nso = [PvSpikeWidth<lmt; SomSpikeWidth<lmt; NonTaggedAllSpikeWidth<lmt];
tags_wso = [PvSpikeWidth>lmt; SomSpikeWidth>lmt; NonTaggedAllSpikeWidth>lmt];

keyboard

% -------------------------------------------------------------------------
function S = smoothallcells(X)

% Smoothing with moving average
sx = size(X);
S = zeros(sx);
for k = 1:sx(1)
    S(k,:) = smooth(X(k,:),'linear',11);
end

% -------------------------------------------------------------------------
function [allpsth5 nbins] = convertpsth(allpsth,template,nbins)

% Input argument check
error(nargchk(2,3,nargin))

% Correlate with template PSTH
allpsth2 = allpsth * mean(template)';

% Optimal bin number
if nargin < 3
    iqr = prctile(allpsth2,75) - prctile(allpsth2,25);
    binsize = 2 * iqr / length(allpsth2) .^ (1 / 3);   % number of levels
    nbins = ceil((max(allpsth2)-min(allpsth2))/binsize);   % Freedman-Diaconis rule
end

% Discretize
allpsth3 = (allpsth2 - min(allpsth2(:))) / (max(allpsth2(:)) - min(allpsth2(:)));   % between 0 and 1
allpsth4 = fix(allpsth3*nbins-1.5) + 2;

% Convert PSTHs to numbers
ep = size(allpsth4,2)-1:-1:0;
tep = 10 .^ ep;
allpsth5 = allpsth4 * tep';

% -------------------------------------------------------------------------
function [msrg srg] = calc_surrogate(allpsth,sss,tags,nbins)

% Calculate surrogate
srg = nan(1,sss);
for k = 1:sss
    
    % Randomize tags
    rr = randperm(length(tags));
    tt = tags(rr);
    
    % Convert PSTHs to numbers
    allpsth5_r = convertpsth(allpsth,allpsth(tt,:),nbins);
    
    % Kullback-Leibler divergence for Pv
    srg(k) = KLdist_discrete(allpsth5_r(tt),allpsth5_r);
end
msrg = median(srg);