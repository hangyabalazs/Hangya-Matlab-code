%% Load & define PSTHs
load('C:\My Dropbox\KepecsLab\_Duda\for balazs\pop_psth_for_PCA\HomeZout.mat')
pv_MatrixPsth1 = pv;
som_MatrixPsth1 = som;
non_tagged_MatrixPsth1 = non_tagged_WS;

load('C:\My Dropbox\KepecsLab\_Duda\for balazs\pop_psth_for_PCA\HomeZin.mat')
pv_MatrixPsth2 = pv;
som_MatrixPsth2 = som;
non_tagged_MatrixPsth2 = non_tagged_WS;

% Load spike width data
load('C:\My Dropbox\KepecsLab\_Duda\for balazs\pop_psth_for_PCA\SpikeWidth.mat')

%% Remove NaNs
[x1 y1] = find(isnan(non_tagged_MatrixPsth1));
[x2 y2] = find(isnan(non_tagged_MatrixPsth2));
x = union(unique(x1),unique(x2));
non_tagged_MatrixPsth1(x,:) = [];
non_tagged_MatrixPsth2(x,:) = [];
NonTaggedSpikeWidth(x) = [];

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
nbins = 1:1200;   % number of levels
for k = nbins
allpsth3 = (allpsth2 - min(allpsth2(:))) / (max(allpsth2(:)) - min(allpsth2(:)));   % between 0 and 1
allpsth4 = fix(allpsth3*k-1.5) + 2;
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
mi(k) = I_X_Yc / H_Xc;

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

er(k) = H_PV / H_all;

kl(k) = KLdist_discrete(allpsth5(1:sizepv(1)),allpsth5);
end

figure;plot(mi,'r')
figure;plot(er,'r')
figure;plot(kl,'r')

%% Spike width

% Discretize
nbins = 1:1200;   % number of levels
for k = nbins
allpsth3 = (allpsth2 - min(allpsth2(:))) / (max(allpsth2(:)) - min(allpsth2(:)));   % between 0 and 1
allpsth4 = fix(allpsth3*k-1.5) + 2;
% allpsth4 = allpsth3 > 0.64;

% Convert PSTHs to numbers
ep = size(allpsth4,2)-1:-1:0;
tep = 10 .^ ep;
allpsth5 = allpsth4 * tep';

% Tag distribution
tags = [PvSpikeWidth<300; SomSpikeWidth<300; NonTaggedSpikeWidth<300];  % knowing spike width

% Mutual information
[I_X_Y H_X H_Y I_X_Yc H_Xc H_Yc] = mutual_information(allpsth5,tags);
mi(k) = I_X_Yc / H_Xc;

H_all = H_Xc;

% Entropy
[I_X_Y H_X H_Y I_X_Yc H_Xc H_Yc] = mutual_information(allpsth5(tags),tags(tags));

H_PV = H_Xc;

er(k) = H_PV / H_all;

kl(k) = KLdist_discrete(allpsth5(tags),allpsth5);
end

figure;plot(mi,'k')
figure;plot(er,'k')
figure;plot(kl,'k')