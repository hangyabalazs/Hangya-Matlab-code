%% Load & define PSTHs
load('C:\My Dropbox\KepecsLab\_Duda\for balazs\pop_psth_for_PCA\HomeZout.mat')
pv_MatrixPsth1 = pv;
som_MatrixPsth1 = som;
non_tagged_MatrixPsth1 = non_tagged_all;

load('C:\My Dropbox\KepecsLab\_Duda\for balazs\pop_psth_for_PCA\HomeZin.mat')
pv_MatrixPsth2 = pv;
som_MatrixPsth2 = som;
non_tagged_MatrixPsth2 = non_tagged_all;

% Load spike width data
load('C:\My Dropbox\KepecsLab\_Duda\for balazs\pop_psth_for_PCA\SpikeWidth.mat')

%% Remove NaNs
[x1 y1] = find(isnan(non_tagged_MatrixPsth1));
[x2 y2] = find(isnan(non_tagged_MatrixPsth2));
x = union(unique(x1),unique(x2));
non_tagged_MatrixPsth1(x,:) = [];
non_tagged_MatrixPsth2(x,:) = [];
NonTaggedAllSpikeWidth(x) = [];

%% Restrict in time
wn1 = -60;
wn2 = 60;
mdl = (size(pv,2) + 1) / 2;
pv_MatrixPsth1 = pv_MatrixPsth1(:,mdl+wn1:mdl+wn2-1);
som_MatrixPsth1 = som_MatrixPsth1(:,mdl+wn1:mdl+wn2-1);
non_tagged_MatrixPsth1 = non_tagged_MatrixPsth1(:,mdl+wn1:mdl+wn2-1);
pv_MatrixPsth2 = pv_MatrixPsth2(:,mdl+wn1:mdl+wn2-1);
som_MatrixPsth2 = som_MatrixPsth2(:,mdl+wn1:mdl+wn2-1);
non_tagged_MatrixPsth2 = non_tagged_MatrixPsth2(:,mdl+wn1:mdl+wn2-1);

%% Smooth
sizepv = size(pv_MatrixPsth1);
spv1 = zeros(sizepv(1),sizepv(2));
for k = 1:sizepv(1)
    spv1(k,:) = smooth(pv_MatrixPsth1(k,:),'linear',11);
end

sizesom = size(som_MatrixPsth1);
ssom1 = zeros(sizesom(1),sizesom(2));
for k = 1:sizesom(1)
    ssom1(k,:) = smooth(som_MatrixPsth1(k,:),'linear',11);
end

sizent = size(non_tagged_MatrixPsth1);
snt1 = zeros(sizent(1),sizent(2));
for k = 1:sizent(1)
    snt1(k,:) = smooth(non_tagged_MatrixPsth1(k,:),'linear',11);
end

spv2 = zeros(sizepv(1),sizepv(2));
for k = 1:sizepv(1)
    spv2(k,:) = smooth(pv_MatrixPsth2(k,:),'linear',11);
end

ssom2 = zeros(sizesom(1),sizesom(2));
for k = 1:sizesom(1)
    ssom2(k,:) = smooth(som_MatrixPsth2(k,:),'linear',11);
end

snt2 = zeros(sizent(1),sizent(2));
for k = 1:sizent(1)
    snt2(k,:) = smooth(non_tagged_MatrixPsth2(k,:),'linear',11);
end

%% Treat segments together
% som_merged = [ssom1 ssom2];
% pv_merged = [spv1 spv2];
% nt_merged = [snt1 snt2];
som_merged = ssom1;
pv_merged = spv1;
nt_merged = snt1;
allpsth = [pv_merged; som_merged; nt_merged];

% Collapse to 10 ms epochs
% wnn = 10;
% rs = reshape(allpsth',wnn,size(allpsth,2)/wnn,size(allpsth,1));
% allpsth2 = squeeze(sum(rs,1))';

% Correlate with mean PV PSTH
allpsth2 = allpsth * mean(pv_merged)';

% Discretize
iqr = prctile(allpsth2,75) - prctile(allpsth2,25);
binsize = 2 * iqr / length(allpsth2) .^ (1 / 3);   % number of levels
nbins = ceil((max(allpsth2)-min(allpsth2))/binsize);   % Freedman-Diaconis rule
nbins1 = nbins;

allpsth3 = (allpsth2 - min(allpsth2(:))) / (max(allpsth2(:)) - min(allpsth2(:)));   % between 0 and 1
allpsth4 = fix(allpsth3*nbins-1.5) + 2;
% allpsth4 = allpsth3 > 0.64;

% Convert PSTHs to numbers
ep = size(allpsth4,2)-1:-1:0;
tep = 10 .^ ep;
allpsth5 = allpsth4 * tep';

% Tag distribution
tags = [zeros(size(pv_merged,1),1); ones(size(som_merged,1)+size(nt_merged,1),1)];  % knowing PV
% tags = [ones(size(pv_merged,1),1); zeros(size(som_merged,1),1); ones(size(nt_merged,1),1)];   % knowing SOM

% Mutual information
[I_X_Y H_X H_Y I_X_Yc H_Xc H_Yc] = mutual_information(allpsth5,tags);
mi = I_X_Yc / H_Xc;

H_all = H_Xc;

% srg = nans(1,50);
% for k = 1:50
%     tags = tags(randperm(length(tags)));
%     [I_X_Y H_X H_Y I_X_Yc H_Xc H_Yc] = mutual_information(allpsth5,tags);
%     srg(k) = I_X_Y / H_X;
% end
% figure
% hist(srg)

% Entropy
[I_X_Y H_X H_Y I_X_Yc H_Xc H_Yc] = mutual_information(allpsth5(1:sizepv(1)),tags(1:sizepv(1)));

H_PV = H_Xc;

er = H_PV / H_all;

kl = KLdist_discrete(allpsth5(1:sizepv(1)),allpsth5);

%% Spike width

L = 180:10:700
next = 1;
for lmt = L

% Tag distribution
tags = [PvSpikeWidth>lmt; SomSpikeWidth>lmt; NonTaggedAllSpikeWidth>lmt];  % knowing spike width

% Correlate with mean PV PSTH
allpsth2 = allpsth * mean(allpsth(tags,:))';

% Discretize
iqr = prctile(allpsth2,75) - prctile(allpsth2,25);
binsize = 2 * iqr / length(allpsth2) .^ (1 / 3);   % number of levels
nbins = ceil((max(allpsth2)-min(allpsth2))/binsize);   % Freedman-Diaconis rule
nbins2 = nbins;

% Discretize
allpsth3 = (allpsth2 - min(allpsth2(:))) / (max(allpsth2(:)) - min(allpsth2(:)));   % between 0 and 1
allpsth4 = fix(allpsth3*nbins-1.5) + 2;
% allpsth4 = allpsth3 > 0.64;

% Convert PSTHs to numbers
ep = size(allpsth4,2)-1:-1:0;
tep = 10 .^ ep;
allpsth5 = allpsth4 * tep';

% Mutual information
[I_X_Y H_X H_Y I_X_Yc H_Xc H_Yc] = mutual_information(allpsth5,tags);
mi(next) = I_X_Yc / H_Xc;

H_all = H_Xc;

% Entropy
[I_X_Y H_X H_Y I_X_Yc H_Xc H_Yc] = mutual_information(allpsth5(tags),tags(tags));

H_PV = H_Xc;

er(next) = H_PV / H_all;

kl(next) = KLdist_discrete(allpsth5(tags),allpsth5);

next = next + 1;
end

figure;
plot(L,kl)

%% significance

% Calculate bootstrap error for PV
B = 200;
kl_bst = nan(1,B);
kl_ns_bst = nan(1,B);
for k = 1:B
    disp(k)
    rrpv = round(rand(1,sizepv(1))*(sizepv(1)-1)+0.5);
    rrsom = round(rand(1,sizesom(1))*(sizesom(1)-1)+0.5);
    rrnt = round(rand(1,sizent(1))*(sizent(1)-1)+0.5);
    
    som_merged = ssom1(rrsom,:);
    pv_merged = spv1(rrpv,:);
    nt_merged = snt1(rrnt,:);
    allpsth = [pv_merged; som_merged; nt_merged];
    tags = [zeros(size(pv_merged,1),1); ones(size(som_merged,1)+size(nt_merged,1),1)];  % knowing PV
    allpsth2 = allpsth * mean(pv_merged)';
    allpsth3 = (allpsth2 - min(allpsth2(:))) / (max(allpsth2(:)) - min(allpsth2(:)));   % between 0 and 1
    allpsth4 = fix(allpsth3*nbins1-1.5) + 2;
    ep = size(allpsth4,2)-1:-1:0;
    tep = 10 .^ ep;
    allpsth5 = allpsth4 * tep';
    kl_bst(k) = KLdist_discrete(allpsth5(1:sizepv(1)),allpsth5);
    
    lmt = 295;
    tags2 = [PvSpikeWidth(rrpv)<lmt; SomSpikeWidth(rrsom)<lmt; NonTaggedSpikeWidth(rrnt)<lmt];  % knowing spike width
    allpsth2 = allpsth * mean(allpsth(tags2,:))';
    allpsth3 = (allpsth2 - min(allpsth2(:))) / (max(allpsth2(:)) - min(allpsth2(:)));   % between 0 and 1
    allpsth4 = fix(allpsth3*nbins2-1.5) + 2;
    ep = size(allpsth4,2)-1:-1:0;
    tep = 10 .^ ep;
    allpsth5 = allpsth4 * tep';
    kl_ns_bst(k) = KLdist_discrete(allpsth5(tags2),allpsth5);
end
kl_err = std(kl_bst);
kl_ns_err = std(kl_ns_bst);

%% plot

figure
bar([0.7414 0.2])
hold on
errorbar([0.7414 0.2],[0.18 0.075])
errorbar([0.7414 0.2],[0.18 0.075],'k+')

%% tests

som_merged = ssom1;
pv_merged = spv1;
nt_merged = snt1;
allpsth = [pv_merged; som_merged; nt_merged];

srg = nan(1,50);
for k = 1:50
    rr = randperm(length(tags));
    tt = tags(rr);
    
    pv_merged = allpsth(tt==0,:);
    others_merged = allpsth(tt==1,:);
    allpsth = [pv_merged; others_merged];
    tags = [zeros(size(pv_merged,1),1); ones(size(som_merged,1)+size(nt_merged,1),1)];  % knowing PV
    allpsth2 = allpsth * mean(pv_merged)';
    allpsth3 = (allpsth2 - min(allpsth2(:))) / (max(allpsth2(:)) - min(allpsth2(:)));   % between 0 and 1
    allpsth4 = fix(allpsth3*nbins1-1.5) + 2;
    ep = size(allpsth4,2)-1:-1:0;
    tep = 10 .^ ep;
    allpsth5 = allpsth4 * tep';
    srg(k) = KLdist_discrete(allpsth5(1:sizepv(1)),allpsth5);
end
figure
hist(srg)

%% test NS

som_merged = ssom1;
pv_merged = spv1;
nt_merged = snt1;
allpsth = [pv_merged; som_merged; nt_merged];
tags2 = [PvSpikeWidth<lmt; SomSpikeWidth<lmt; NonTaggedSpikeWidth<lmt];  % knowing spike width

srg = nan(1,50);
for k = 1:50
    rr = randperm(length(tags2));
    tt = tags2(rr);
    
    tags2 = tt;  % knowing spike width
    allpsth2 = allpsth * mean(allpsth(tags2,:))';
    allpsth3 = (allpsth2 - min(allpsth2(:))) / (max(allpsth2(:)) - min(allpsth2(:)));   % between 0 and 1
    allpsth4 = fix(allpsth3*nbins2-1.5) + 2;
    ep = size(allpsth4,2)-1:-1:0;
    tep = 10 .^ ep;
    allpsth5 = allpsth4 * tep';
    srg(k) = KLdist_discrete(allpsth5(tags2),allpsth5);
end
figure
hist(srg)

%% check for effect of sample size on KL

next = 1;
for k = 1:1177
    tags = [false(1,k) true(1,1177-k)];
    tags = tags(randperm(1177));
    % Correlate with mean PV PSTH
allpsth2 = allpsth * mean(allpsth)';

% Discretize
iqr = prctile(allpsth2,75) - prctile(allpsth2,25);
binsize = 2 * iqr / length(allpsth2) .^ (1 / 3);   % number of levels
nbins = ceil((max(allpsth2)-min(allpsth2))/binsize);   % Freedman-Diaconis rule
nbins2 = nbins;
nbins = 10;

% Discretize
allpsth3 = (allpsth2 - min(allpsth2(:))) / (max(allpsth2(:)) - min(allpsth2(:)));   % between 0 and 1
allpsth4 = fix(allpsth3*nbins-1.5) + 2;
% allpsth4 = allpsth3 > 0.64;

% Convert PSTHs to numbers
ep = size(allpsth4,2)-1:-1:0;
tep = 10 .^ ep;
allpsth5 = allpsth4 * tep';

kl(next) = KLdist_discrete(allpsth5(tags),allpsth5);

next = next + 1;
end