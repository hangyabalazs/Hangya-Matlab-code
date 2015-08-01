function somentropy_fcn

% Load & define PSTHs
global DATADIR
load([DATADIR 'SOM_Sachin\PSTH\HomeZout.mat'])
pv_MatrixPsth1 = pv;
som_MatrixPsth1 = som;
non_tagged_MatrixPsth1 = non_tagged_all;

load([DATADIR 'SOM_Sachin\PSTH\HomeZin.mat'])
pv_MatrixPsth2 = pv;
som_MatrixPsth2 = som;
non_tagged_MatrixPsth2 = non_tagged_all;

% Load spike width data
load([DATADIR 'SOM_Sachin\PSTH\SpikeWidth.mat'])

% Remove NaNs
[x1 y1] = find(isnan(non_tagged_MatrixPsth1));
[x2 y2] = find(isnan(non_tagged_MatrixPsth2));
x = union(unique(x1),unique(x2));
non_tagged_MatrixPsth1(x,:) = [];
non_tagged_MatrixPsth2(x,:) = [];
NonTaggedAllSpikeWidth(x) = [];

% Restrict in time
wn1 = -60;
wn2 = 60;
mdl = (size(pv,2) + 1) / 2;
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
kl_pv = KLdist_discrete(allpsth5_pv(tags_pv),allpsth5_pv);

% Narrow spiking cells
L = 180:10:380;
kl_ns = nan(size(L));
next = 1;
for lmt = L
    
    % Tag distribution
    tags_ns = [PvSpikeWidth<lmt; SomSpikeWidth<lmt; NonTaggedAllSpikeWidth<lmt];  % knowing narrow spiking
    
    % Correlate with mean PV PSTH
    [allpsth5_ns nbins_ns] = convertpsth(allpsth,allpsth(tags_ns,:));
    
    % Kullback-Leibler divergence
    kl_ns(next) = KLdist_discrete(allpsth5_ns(tags_ns),allpsth5_ns);
    next = next + 1;
end

% Wide spiking cells
kl_ws = nan(size(L));
next = 1;
for lmt = L
    
    % Tag distribution
    tags_ws = [PvSpikeWidth>lmt; SomSpikeWidth>lmt; NonTaggedAllSpikeWidth>lmt];  % knowing wide spiking
    
    % Correlate with mean PV PSTH
    [allpsth5_ws nbins_ws] = convertpsth(allpsth,allpsth(tags_ws,:));
    
    % Kullback-Leibler divergence
    kl_ws(next) = KLdist_discrete(allpsth5_ws(tags_ws),allpsth5_ws);
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
tags_nso = [PvSpikeWidth<lmt; SomSpikeWidth<lmt; NonTaggedAllSpikeWidth<lmt];
tags_wso = [PvSpikeWidth>lmt; SomSpikeWidth>lmt; NonTaggedAllSpikeWidth>lmt];

% Calculate bootstrap error
B = 200;
kl_pv_bst = nan(1,B);
kl_ns_bst = nan(1,B);
for k = 1:B
    
    % Resample PV, SOM and non-tagged populations
    rrpv = round(rand(1,sizepv(1))*(sizepv(1)-1)+0.5);
    rrsom = round(rand(1,sizesom(1))*(sizesom(1)-1)+0.5);
    rrnt = round(rand(1,sizent(1))*(sizent(1)-1)+0.5);
    
    pv_merged_r = spv1(rrpv,:);
    som_merged_r = ssom1(rrsom,:);
    nt_merged_r = snt1(rrnt,:);
    allpsth_r = [pv_merged_r; som_merged_r; nt_merged_r];
    
    % Convert PSTHs to numbers
    allpsth5_pv_r = convertpsth(allpsth_r,pv_merged_r,nbins_pv);
    
    % Kullback-Leibler divergence for PV
    kl_pv_bst(k) = KLdist_discrete(allpsth5_pv_r(tags_pv),allpsth5_pv_r);
    
    % Convert PSTHs to numbers
    tags_ns_r = [PvSpikeWidth(rrpv)<lmt; SomSpikeWidth(rrsom)<lmt; NonTaggedAllSpikeWidth(rrnt)<lmt];  % knowing spike width
    allpsth5_pv_r = convertpsth(allpsth_r,allpsth_r(tags_ns_r),nbins_ns);
    
    % Kullback-Leibler divergence for NS
    kl_ns_bst(k) = KLdist_discrete(allpsth5_pv_r(tags_ns_r),allpsth5_pv_r);
end
kl_pv_err = std(kl_pv_bst);
kl_ns_err = std(kl_ns_bst);

% Plot
figure
bar([kl_pv kl_nso])
hold on
errorbar([kl_pv kl_nso],[kl_pv_err kl_ns_err],'k+')

% Test whether KL for PV is significant
srg = nan(1,50);
for k = 1:50
    
    % Randomize tags
    rr = randperm(length(tags_pv));
    tt = tags_pv(rr);
    pv_merged_r = allpsth(tt,:);
    others_merged_r = allpsth(~tt,:);
    allpsth_r = [pv_merged_r; others_merged_r];
    
    % Convert PSTHs to numbers
    allpsth5_pv_r = convertpsth(allpsth_r,pv_merged_r,nbins_pv);
    
    % Kullback-Leibler divergence for Pv
    srg(k) = KLdist_discrete(allpsth5_pv_r(tags_pv),allpsth5_pv_r);
end
figure
hist(srg)
p_pv = sum(srg>kl_pv) / length(srg);
disp(p_pv)

% Test whether KL for NS is significant
srg = nan(1,50);
for k = 1:50
    
    % Randomize tags
    rr = randperm(length(tags_nso));
    tt = tags_nso(rr);
    
    % Convert PSTHs to numbers
    allpsth5_ns_r = convertpsth(allpsth,allpsth(tt,:),nbins_ns);
    
    % Kullback-Leibler divergence for Pv
    srg(k) = KLdist_discrete(allpsth5_ns_r(tt),allpsth5_ns_r);
end
figure
hist(srg)
p_ns = sum(srg>kl_nso) / length(srg);
disp(p_ns)

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