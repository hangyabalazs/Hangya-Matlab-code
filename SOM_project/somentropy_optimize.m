function somentropy_optimize
%SOMENTROPY_OPTIMIZE   Kullbach-Leibler divergence for PSTH in different groups.
%   SOMENTROPY_OPTIMIZE linearizes PSTH distributions by taking the
%   correlation of each PSTH with the corresponding group mean (e.g. mean
%   PV PSTH). PSTHs are restricted to -600 - +600 ms windows around an event
%   (e.g. HomeZoneOut by default) and smoothed with a moving average.
%
%   Then, PSTH correlation values are discretized using the
%   Freedman-Diaconis rule and Kullbach-Leibler distance between PV PSTH
%   and all PSTH correlation values are calculated.
%
%   Similarly, KL distance for narrow spiking and wide spiking cell groups
%   are calculated. Spike width threshold is optimized separately for NS
%   and WS between 210 and 350 us to get maximal KL value.
%
%   Bootstrap errors are computed by generating a 200 bootrsrap sample,
%   resampling the PSTH distributions for the various groups.
%
%   Significance of KL distances as well as significant differences among
%   them are assesed by randomly reassigning the group tags.
%
%   See also KLDIST_DISCRETE and KLDIST.

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
load([DATADIR 'SOM_Sachin\PSTH\Indexes2Exclude.mat'])

% Exlude
% pv_MatrixPsth1(ind_excl_pv,:) = [];
% pv_MatrixPsth2(ind_excl_pv,:) = [];
% PvSpikeWidth(ind_excl_pv) = [];
% non_tagged_MatrixPsth1(ind_excl_all,:) = [];
% non_tagged_MatrixPsth2(ind_excl_all,:) = [];
% NonTaggedAllSpikeWidth(ind_excl_all) = [];

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
L = 210:10:350;
kl_ns = nan(size(L));
next = 1;
for lmt = L
    disp(lmt)
    
    % Tag distribution
    tags_ns = [PvSpikeWidth<lmt; SomSpikeWidth<lmt; NonTaggedAllSpikeWidth<lmt];  % knowing narrow spiking
    
    % Correlate with mean PV PSTH
    [allpsth5_ns nbins_ns] = convertpsth(allpsth,allpsth(tags_ns,:),nbins_pv);   % note: nbins_pv is used for NS and WS as well
    
    % Kullback-Leibler divergence
    kl_ns(next) = KLdist_discrete(allpsth5_ns(tags_ns),allpsth5_ns);
    
    next = next + 1;
end

% Wide spiking cells
kl_ws = nan(size(L));
next = 1;
for lmt = L
    disp(lmt)
    
    % Tag distribution
    tags_ws = [PvSpikeWidth>lmt; SomSpikeWidth>lmt; NonTaggedAllSpikeWidth>lmt];  % knowing wide spiking
    
    % Correlate with mean PV PSTH
    [allpsth5_ws nbins_ws] = convertpsth(allpsth,allpsth(tags_ws,:),nbins_pv);
    
    % Kullback-Leibler divergence
    kl_ws(next) = KLdist_discrete(allpsth5_ws(tags_ws),allpsth5_ws);
    
    next = next + 1;
end

% Optimize cut for spike width
% sumkl = kl_ns + kl_ws;
figure;
plot(L,kl_ns,'Color',[0.8 0.2 0.5])
hold on
plot(L,kl_ws,'Color',[0.2 0.8 0.5])
inx = kl_ns == max(kl_ns);
lmt_ns = L(inx);
kl_nso = kl_ns(inx);
tags_nso = [PvSpikeWidth<lmt_ns; SomSpikeWidth<lmt_ns; NonTaggedAllSpikeWidth<lmt_ns];
inx = kl_ws == max(kl_ws);
lmt_ws = L(inx);
kl_wso = kl_ws(inx);
tags_wso = [PvSpikeWidth>lmt_ws; SomSpikeWidth>lmt_ws; NonTaggedAllSpikeWidth>lmt_ws];

% Calculate bootstrap error
B = 200;
kl_pv_bst = nan(1,B);
kl_ns_bst = nan(1,B);
kl_ws_bst = nan(1,B);
for k = 1:B
    disp(k)
    
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
    tags_ns_r = [PvSpikeWidth(rrpv)<lmt_ns; SomSpikeWidth(rrsom)<lmt_ns; NonTaggedAllSpikeWidth(rrnt)<lmt_ns];  % knowing spike width
    allpsth5_ns_r = convertpsth(allpsth_r,allpsth_r(tags_ns_r),nbins_pv);
    
    % Kullback-Leibler divergence for NS
    kl_ns_bst(k) = KLdist_discrete(allpsth5_ns_r(tags_ns_r),allpsth5_ns_r);
    
    % Convert PSTHs to numbers
    tags_ws_r = [PvSpikeWidth(rrpv)>lmt_ws; SomSpikeWidth(rrsom)>lmt_ws; NonTaggedAllSpikeWidth(rrnt)>lmt_ws];  % knowing spike width
    allpsth5_ws_r = convertpsth(allpsth_r,allpsth_r(tags_ws_r),nbins_pv);
    
    % Kullback-Leibler divergence for NS
    kl_ws_bst(k) = KLdist_discrete(allpsth5_ws_r(tags_ws_r),allpsth5_ws_r);
end
kl_pv_err = std(kl_pv_bst);
kl_ns_err = std(kl_ns_bst);
kl_ws_err = std(kl_ws_bst);

% Plot
figure
bar([kl_pv kl_nso kl_wso])
hold on
errorbar([kl_pv kl_nso kl_wso],[kl_pv_err kl_ns_err kl_ws_err],'k+')

% Test whether KL for PV is significant
sss = 100;
[msrg srg] = calc_surrogate(allpsth,sss,tags_pv,nbins_pv);
figure
hist(srg)
p_pv = sum(srg>kl_pv) / length(srg);
disp(p_pv)

% Test whether KL for NS is significant
[msrg srg] = calc_surrogate(allpsth,sss,tags_nso,nbins_pv);
figure
hist(srg)
p_ns = sum(srg>kl_nso) / length(srg);
disp(p_ns)

% Test whether KL for WS is significant
[msrg srg] = calc_surrogate(allpsth,sss,tags_wso,nbins_pv);
figure
hist(srg)
p_ws = sum(srg>kl_wso) / length(srg);
disp(p_ws)

% Compare PV ans NS
sss2 = 1000;
srg_pv_ns = nan(1,sss2);
srg_pv_ws = nan(1,sss2);
srg_ns_ws = nan(1,sss2);
srg_pv = nan(1,sss2);
srg_ns = nan(1,sss2);
srg_ws = nan(1,sss2);
for k = 1:sss2
    
    % Randomize tags - PV
    rr = randperm(length(tags_pv));
    tt = tags_pv(rr);
    
    % Convert PSTHs to numbers
    allpsth5_r = convertpsth(allpsth,allpsth(tt,:),nbins_pv);
    
    % Kullback-Leibler divergence for PV
    srg_pv(k) = KLdist_discrete(allpsth5_r(tt),allpsth5_r);
    
    % Randomize tags - NS
    rr = randperm(length(tags_nso));
    tt = tags_nso(rr);
    
    % Convert PSTHs to numbers
    allpsth5_r = convertpsth(allpsth,allpsth(tt,:),nbins_pv);
    
    % Kullback-Leibler divergence for NS
    srg_ns(k) = KLdist_discrete(allpsth5_r(tt),allpsth5_r);
    
    % Randomize tags - WS
    rr = randperm(length(tags_wso));
    tt = tags_wso(rr);
    
    % Convert PSTHs to numbers
    allpsth5_r = convertpsth(allpsth,allpsth(tt,:),nbins_pv);
    
    % Kullback-Leibler divergence for NS
    srg_ws(k) = KLdist_discrete(allpsth5_r(tt),allpsth5_r);
    
    % Differences
    srg_pv_ns(k) = srg_pv(k) - srg_ns(k);
    srg_pv_ws(k) = srg_pv(k) - srg_ws(k);
    srg_ns_ws(k) = srg_ns(k) - srg_ws(k);
end
figure
hist(srg_pv_ns)
p_pv_ns = sum(srg_pv_ns>(kl_pv-kl_nso)) / length(srg_pv_ns);
disp(p_pv_ns)       % p = 0.035; 0.003 with wn = 30ms; 0.008 with excusion (wn = 60 ms)
p_pv_ws = sum(srg_pv_ws>(kl_pv-kl_wso)) / length(srg_pv_ws);
disp(p_pv_ws)       % p = 0.002
p_ns_ws = sum(srg_ns_ws>(kl_nso-kl_wso)) / length(srg_ns_ws);
disp(p_ns_ws)       % p = 0.002
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