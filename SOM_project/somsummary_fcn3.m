function somsummary_fcn3
%SOMSUMMARY_FCN3   Linear discriminant analysis.
%   SOMSUMMARY_FCN3 calculates LDA (linear discriminant analysis) on PSTHs
%   of identified interneurons to show separation of subtypes (as summary
%   result). LDA is calculated on three groups - identified PV and SOM
%   cells and concurrently recorded wide-spiking cells. Seven dimensions
%   are used, correlation with mean PV and SOM PSTHs aligned to HomeZoneIn 
%   and HomeZoneOut events (4 dimensions), supression index for HomeZoneIn
%   and supression time for PSTHs aligned to both events (2 dimensions).
%   Supression index is difference of firing between second and middle
%   portion of the PSTH. Supression time is time of the PSTH with blow 25%
%   firing. Extremes are estimated as 1 and 99 percentiles for each group.
%
%   LDA returns the plane containing the centroids of the three included
%   groups. Since this plane is unique, no additional steps is taken. (Note
%   that for more groups, LDA finds dimensions of the subspace containing
%   the centroids according to the variance of the groups. This does not
%   happen if the number of groups is only 3.) All four groups as well as
%   the group of all unidentified cells is projected onto this plane. An
%   additional PCA is carried out on the returning dimensions to remove
%   correlations between the axes.
%
%   See also SOMENTROPY_OPTIMIZE.

% Load HomeZoneOut data
global DATADIR
load([DATADIR 'SOM_Sachin\PSTH\HomeZout.mat'])
pv_MatrixPsth1 = pv_psth;
som_MatrixPsth1 = som_psth;
non_tagged_MatrixPsth1 = nt_psth;
tetrode_pairs_WS1 = WS_tetrode_pairs_psth;

% Load HomeZoneIn data
load([DATADIR 'SOM_Sachin\PSTH\HomeZin.mat'])
pv_MatrixPsth2 = pv_psth;
som_MatrixPsth2 = som_psth;
non_tagged_MatrixPsth2 = nt_psth;
tetrode_pairs_WS2 = WS_tetrode_pairs_psth;

% Remove NaNs
[x1 y1] = find(isnan(non_tagged_MatrixPsth1));
[x2 y2] = find(isnan(non_tagged_MatrixPsth2));
x = union(unique(x1),unique(x2));
non_tagged_MatrixPsth1(x,:) = [];
non_tagged_MatrixPsth2(x,:) = [];

% Restrict in time
pv_MatrixPsth1 = pv_MatrixPsth1(:,7:207);
som_MatrixPsth1 = som_MatrixPsth1(:,7:207);
non_tagged_MatrixPsth1 = non_tagged_MatrixPsth1(:,7:207);
tetrode_pairs_WS1 = tetrode_pairs_WS1(:,7:207);
pv_MatrixPsth2 = pv_MatrixPsth2(:,7:207);
som_MatrixPsth2 = som_MatrixPsth2(:,7:207);
non_tagged_MatrixPsth2 = non_tagged_MatrixPsth2(:,7:207);
tetrode_pairs_WS2 = tetrode_pairs_WS2(:,7:207);

% Smooth
sizepv = size(pv_MatrixPsth1);
sizesom = size(som_MatrixPsth1);
spv1 = smoothallcells(pv_MatrixPsth1);   % PV, HomeZoneOut
ssom1 = smoothallcells(som_MatrixPsth1);   % SOM, HomeZoneOut
snt1 = smoothallcells(non_tagged_MatrixPsth1);   % non-tagged, HomeZoneOut
stp1 = smoothallcells(tetrode_pairs_WS1);   % tetrode pairs, HomeZoneOut
spv2 = smoothallcells(pv_MatrixPsth2);   % PV, HomeZoneIn
ssom2 = smoothallcells(som_MatrixPsth2);   % SOM, HomeZoneIn
snt2 = smoothallcells(non_tagged_MatrixPsth2);   % non-tagged, HomeZoneIn
stp2 = smoothallcells(tetrode_pairs_WS2);   % tetrode pairs, HomeZoneIn

% Supression index (for HomeZoneIn); difference of firing in second half
% and mid portion of the segment
inx2 = 111:200;    % second half of the segment
inx1 = 80:110;     % middle portion of the segment
si_pv2 = supression_index(spv2,inx1,inx2);    % PV
si_som2 = supression_index(ssom2,inx1,inx2);  % SOM
si_stp2 = supression_index(stp2,inx1,inx2);   % tetrode pairs
si_snt2 = supression_index(snt2,inx1,inx2);   % non-tagged

% Correlation with specific PSTHs
template = mean(spv1);   % correlation with mean PV, HomeZoneOut
spvd1 = tcorr(spv1,template);    % PV, HomeZoneOut
ssomd1 = tcorr(ssom1,template);  % SOM, HomeZoneOut
stpd1 = tcorr(stp1,template);    % tetrode pairs, HomeZoneOut
sntd1 = tcorr(snt1,template);    % non-tagged, HomeZoneOut

template2 = mean(ssom1);    % correlation with mean SOM, HomeZoneOut
spve1 = tcorr(spv1,template2);   % PV, HomeZoneOut
ssome1 = tcorr(ssom1,template2); % SOM, HomeZoneOut
stpe1 = tcorr(stp1,template2);   % tetrode pairs, HomeZoneOut
snte1 = tcorr(snt1,template2);   % non-tagged, HomeZoneOut

template = mean(spv2);   % correlation with mean PV, HomeZoneIn
spvd2 = tcorr(spv2,template);    % PV, HomeZoneIn
ssomd2 = tcorr(ssom2,template);  % SOM, HomeZoneIn
stpd2 = tcorr(stp2,template);    % tetrode pairs, HomeZoneIn
sntd2 = tcorr(snt2,template);    % non-tagged, HomeZoneIn

template2 = mean(ssom2);   % correlation with mean SOM, HomeZoneIn
spve2 = tcorr(spv2,template2);   % PV, HomeZoneIn
ssome2 = tcorr(ssom2,template2); % SOM, HomeZoneIn
stpe2 = tcorr(stp2,template2);   % tetrode pairs, HomeZoneIn
snte2 = tcorr(snt2,template2);   % non-tagged, HomeZoneIn

% Supression time
mlt = 0.25;   % 0.25: maybe best; 0.3 also good
p1 = 1;
p2 = 99;
los_pv1 = supression_time(spv1,mlt,p1,p2);      % PV, HomeZoneOut
los_som1 = supression_time(ssom1,mlt,p1,p2);    % SOM, HomeZoneOut
los_stp1 = supression_time(stp1,mlt,p1,p2);     % tetrode pairs, HomeZoneOut
los_nt1 = supression_time(snt1,mlt,p1,p2);      % non-tagged, HomeZoneOut

los_pv2 = supression_time(spv2,mlt,p1,p2);      % PV, HomeZoneIn
los_som2 = supression_time(ssom2,mlt,p1,p2);    % SOM, HomeZoneIn
los_stp2 = supression_time(stp2,mlt,p1,p2);     % tetrode pairs, HomeZoneIn
los_nt2 = supression_time(snt2,mlt,p1,p2);      % non-tagged, HomeZoneIn

% LDA - PV-pairs and SOM-pairs included as one group
v1 = [spvd1 ssomd1 stpd1]';  % correlation with mean PV, HomeZoneOut
meanv1 = mean(v1);
stdv1 = std(v1);
v1 = (v1 - meanv1) / stdv1;
v1b = [spvd2 ssomd2 stpd2]'; % correlation with mean PV, HomeZoneIn
meanv1b = mean(v1b);
stdv1b = std(v1b);
v1b = (v1b - meanv1b) / stdv1b;
v2 = [spve1 ssome1 stpe1]';  % correlation with mean SOM, HomeZoneOut
meanv2 = mean(v2);
stdv2 = std(v2);
v2 = (v2 - meanv2) / stdv2;
v2b = [spve2 ssome2 stpe2]';  % correlation with mean SOM, HomeZoneIn
meanv2b = mean(v2b);
stdv2b = std(v2b);
v2b = (v2b - meanv2b) / stdv2b;
v3 = [si_pv2 si_som2 si_stp2]';   % supression index, HomeZoneIn
meanv3 = mean(v3);
stdv3 = std(v3);
v3 = (v3 - meanv3) / stdv3;
v4 = [los_pv1 los_som1 los_stp1]';   % supression time, HomeZoneOut
meanv4 = mean(v4);
stdv4 = std(v4);
v4 = (v4 - meanv4) / stdv4;
v5 = [los_pv2 los_som2 los_stp2]';   % supression time, HomeZoneIn
meanv5 = mean(v5);
stdv5 = std(v5);
v5 = (v5 - meanv5) / stdv5;
vars = [v1 v2 v1b v2b...   % corr. w PV, corr. w SOM.
    v3...    % supression index (HZIn)
    v4 v5];  % supression time (HZOut), supression time (HZIn)
celltyps = [zeros(sizepv(1),1); ones(sizesom(1),1); ones(size(stp1,1),1)*2;];

W = lda2(vars,celltyps);   % linear discriminant analysis
u1 = [spvd1 ssomd1 stpd1 sntd1]';
u1 = (u1 - meanv1) / stdv1;   % standardize using the original group means and SDs
u1b = [spvd2 ssomd2 stpd2 sntd2]';
u1b = (u1b - meanv1b) / stdv1b;
u2 = [spve1 ssome1 stpe1 snte1]';
u2 = (u2 - meanv2) / stdv2;
u2b = [spve2 ssome2 stpe2 snte2]';
u2b = (u2b - meanv2b) / stdv2b;
u3 = [si_pv2 si_som2 si_stp2 si_snt2]';
u3 = (u3 - meanv3) / stdv3;
u4 = [los_pv1 los_som1 los_stp1 los_nt1]';
u4 = (u4 - meanv4) / stdv4;
u5 = [los_pv2 los_som2 los_stp2 los_nt2]';
u5 = (u5 - meanv5) / stdv5;
vars2 = [u1 u2 u1b u2b...   % corr. w PV, corr. w SOM.
    u3...    % supression index (HZIn)
    u4 u5];  % supression time (HZOut), supression time (HZIn)

scores = [ones(size(vars2,1),1) vars2] * W';    % project five groups (including all non-tagged) to the plane found
score_1 = scores(:,2);
score_2 = scores(:,3);

celltypes = {1:sizepv(1); sizepv(1)+1:sizepv(1)+sizesom(1); ...
    sizepv(1)+sizesom(1)+1:sizepv(1)+sizesom(1)+size(stp1,1); ...
    sizepv(1)+sizesom(1)+size(stp1,1)+1:size(vars2,1)};

figure    % plot
plot(score_1(celltypes{4}),score_2(celltypes{4}),'.','MarkerSize',15,'Color',[0.7 0.7 0.7]);
hold on
plot(score_1(celltypes{3}),score_2(celltypes{3}),'k.','MarkerSize',20);
plot(score_1(celltypes{1}),score_2(celltypes{1}),'r.','MarkerSize',20);
plot(score_1(celltypes{2}),score_2(celltypes{2}),'b.','MarkerSize',20);

[coeff PCAscores] = princomp([score_1 score_2]);    % remove correlations with PCA
PC1 = PCAscores(:,1);
PC2 = PCAscores(:,2);

figure   % plot
plot(PC1(celltypes{4}),PC2(celltypes{4}),'o','MarkerSize',5,'Color',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7]);
hold on
plot(PC1(celltypes{3}),PC2(celltypes{3}),'ko','MarkerSize',12,'MarkerFaceColor','k');
plot(PC1(celltypes{1}),PC2(celltypes{1}),'ro','MarkerSize',12,'MarkerFaceColor','r');
plot(PC1(celltypes{2}),PC2(celltypes{2}),'bo','MarkerSize',12,'MarkerFaceColor','b');

% -------------------------------------------------------------------------
function S = smoothallcells(X)

% Smoothing with moving average
sx = size(X);
S = zeros(sx);
for k = 1:sx(1)
    S(k,:) = smooth(X(k,:),'linear',11);
end

% -------------------------------------------------------------------------
function S = supression_index(X,inx1,inx2)

% Difference of firing in second half and mid portion of the segment
S = mean(X(:,inx2),2)' - mean(X(:,inx1),2)';

% -------------------------------------------------------------------------
function C = tcorr(X,template)

% Correlation with specific PSTHs
C = template * X';

% -------------------------------------------------------------------------
function ST = supression_time(X,mlt,p1,p2)

% Time of the PSTH with blow 25% firing 
bs = (prctile(X(:),p2) - prctile(X(:),p1)) * mlt + prctile(X(:),p1);   % extremes are estimated as 1 and 99 percentiles for each group
ST = sum(X<bs,2)';