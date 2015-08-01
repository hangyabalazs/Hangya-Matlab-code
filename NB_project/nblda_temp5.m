function nblda
%NBLDA   Linear discriminant analysis.
%   NBLDA calculates LDA (linear discriminant analysis) on PSTHs
%   of identified basal forebrain cells to show separation of subtypes. LDA is 
%   calculated on two groups - identified ChAT and PV
%   cells. Two dimensions are used, correlation with mean ChAT and PV PSTHs
%   aligned to LeftPortIn and z-scored.
%
%   LDA returns the line containing the centroids of the two included
%   groups. Both groups as well as
%   the group of all unidentified cells are projected onto this line. LDA
%   scores are plotted against baseline firing rate.
%
%   See also SOMSUMMARY_FCN2.

% Edit log: BH 8/31/12

% Load data
global DATAPATH
load([DATAPATH 'NB\responsesorter3\allPSTH.mat'])

% PSTH statistics
activation_peak = [allstats.activation_peak];   % peak time of activation
activation_start = [allstats.activation_start];   % activation onset time
activation_end = [allstats.activation_end];   % activation offset time
maxvalue = [allstats.maxvalue];   % maximal firing rate
Wpa = [allstats.Wpa];   % Mann-Whitney test for significant activation
inhibition_peak = [allstats.inhibition_peak];   % peak time of inhibition
inhibition_start = [allstats.inhibition_start];   % inhibition onset time
inhibition_end = [allstats.inhibition_end];   % inhibition offset time
minvalue = [allstats.minvalue];   % minimla firing rate
Wpi = [allstats.Wpi];   % Mann-Whitney test for significant inhibition
baseline = [allstats.baseline];   % baseline firing rate

% Groups of activated and inhibited cells
inx_act = Wpa < 0.01;   % significant activation
inx_inh = Wpi < 0.01;   % significant inhibition
activated = find(inx_act&~inx_inh);    % indices for activated cells
inhibited = find(inx_inh&~inx_act);    % indices for inhibited cells
ai = find(inx_act&inx_inh);
inx = activation_peak(ai) < inhibition_peak(ai);
activated_inhibited = sort(ai(inx));
inhibited_activated = ai(~inx);
activated = [activated'; activated_inhibited'];   % categorize based on first significant effect
inhibited = [inhibited'; inhibited_activated'];

% Identified NB cells
selstr = ['"PV+"==1&"ID_PC">20&"Lr_PC"<0.15&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
PV = selectcell(selstr);    % cell IDs for PV cells
selstr = ['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
ChAT = selectcell(selstr);   % cell IDs for ChAT cells

[nms inxa PVinx] = intersect(PV,tags);   % indices for PV cells
[nms inxa ChATinx] = intersect(ChAT,tags);   % indices for ChAT cells
PVactinx = intersect(PVinx,activated);   % indices for activated PV cells
ChATactinx = intersect(ChATinx,activated);   % indices for activated ChAT cells
NTactinx = setdiff(activated,union(ChATactinx,PVactinx));  % indices for unidentified cells
numChAT = length(ChATactinx);   % number of ChAT cells
numPV = length(PVactinx);   % number of PV cells
numNT = length(NTactinx);   % number of unidentified cells
ChAT_PSTH = allpsth(ChATactinx,:);   % PSTHs of ChAT cells
PV_PSTH = allpsth(PVactinx,:);    % PSTHs of PV cells

% Correlation with specific PSTHs
template = mean(ChAT_PSTH);   % mean PSTH of ChAT cells
ChATcorr = tcorr(allpsth,template);   % correlation with mean ChAT PSTH
template = mean(PV_PSTH);   % mean PSTH of ChAT cells
PVcorr = tcorr(allpsth,template);   % correlation with mean PV PSTH

% Variables
v01 = activation_peak;   % peak time of effect
v02 = log10(maxvalue./baseline);   % relative firing rate change on log scale
v03 = activation_end - activation_start;   % activation duration
v04 = baseline;   % base firing rate
v05 = ChATcorr;   % correlation with mean ChAT PSTH
v06 = PVcorr;   % correlation with mean PV PSTH

% Restrict to ChAT and PV groups
v1 = [v01(ChATactinx) v01(PVactinx)]';
v2 = [v02(ChATactinx) v02(PVactinx)]';
v3 = [v03(ChATactinx) v03(PVactinx)]';
v4 = [v04(ChATactinx) v04(PVactinx)]';
v5 = [v05(ChATactinx) v05(PVactinx)]';
v6 = [v06(ChATactinx) v06(PVactinx)]';

% Standardize variables
meanv1 = mean(v1);
stdv1 = std(v1);
v1 = (v1 - meanv1) / stdv1;
meanv2 = mean(v2);
stdv2 = std(v2);
v2 = (v2 - meanv2) / stdv2;
meanv3 = mean(v3);
stdv3 = std(v3);
v3 = (v3 - meanv3) / stdv3;
meanv4 = mean(v4);
stdv4 = std(v4);
v4 = (v4 - meanv4) / stdv4;
meanv5 = mean(v5);
stdv5 = std(v5);
v5 = (v5 - meanv5) / stdv5;
meanv6 = mean(v6);
stdv6 = std(v6);
v6 = (v6 - meanv6) / stdv6;

% LDA
vars = [v5 v6];
celltyps = [zeros(numChAT,1); ones(numPV,1)];
W = lda2(vars,celltyps);   % linear discriminant analysis

% Add non-identified cells
u1 = [v01(ChATactinx) v01(PVactinx) v01(NTactinx)]';
u2 = [v02(ChATactinx) v02(PVactinx) v02(NTactinx)]';
u3 = [v03(ChATactinx) v03(PVactinx) v03(NTactinx)]';
u4 = [v04(ChATactinx) v04(PVactinx) v04(NTactinx)]';
u5 = [v05(ChATactinx) v05(PVactinx) v05(NTactinx)]';
u6 = [v06(ChATactinx) v06(PVactinx) v06(NTactinx)]';

% 'Standardize' new variables
u1 = (u1 - meanv1) / stdv1;   % standardize using the original group means and SDs
u2 = (u2 - meanv2) / stdv2;
u3 = (u3 - meanv3) / stdv3;
u4 = (u4 - meanv4) / stdv4;
u5 = (u5 - meanv5) / stdv5;
u6 = (u6 - meanv6) / stdv6;

% Plot LDA 
vars2 = [u5 u6];
scores = [ones(size(vars2,1),1) vars2] * W';    % project all groups (including all non-tagged) to the subspace found
score_1 = scores(:,1);
score_2 = scores(:,2);
celltypes = {1:numChAT; numChAT+1:numChAT+numPV; numChAT+numPV+1:size(vars2,1)};

figure    % plot
plot(score_1(celltypes{3}),score_2(celltypes{3}),'.','MarkerSize',15,'Color',[0.7 0.7 0.7]);
hold on
plot(score_1(celltypes{1}),score_2(celltypes{1}),'.','Color',[0 0.8 0],'MarkerSize',20);
plot(score_1(celltypes{2}),score_2(celltypes{2}),'.','Color',[0.8 0 0],'MarkerSize',20);

% PCA
[coeff PCAscores] = princomp([score_1 score_2]);    % remove correlations with PCA
PC1 = PCAscores(:,1);   % for two groups, PC1 is aligned with the line found by LDA
PC2 = PCAscores(:,2);

figure   % plot
plot(PC1(celltypes{3}),PC2(celltypes{3}),'o','MarkerSize',5,'Color',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7]);
hold on
plot(PC1(celltypes{1}),PC2(celltypes{1}),'o','MarkerSize',12,'MarkerFaceColor',[0 0.8 0],'MarkerEdgeColor',[0 0.8 0]);
plot(PC1(celltypes{2}),PC2(celltypes{2}),'o','MarkerSize',12,'MarkerFaceColor',[0.8 0 0],'MarkerEdgeColor',[0.8 0 0]);

% Plot LDA vs baseline firing rate
% u4 = [v04(ChATactinx) v04(PVactinx) v04(NTactinx)]';  % for non-normalized firing rate
figure   % plot
plot(PC1(celltypes{3}),u4(celltypes{3}),'o','MarkerSize',5,'Color',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7]);
hold on
plot(PC1(celltypes{1}),u4(celltypes{1}),'o','MarkerSize',12,'MarkerFaceColor',[0 0.8 0],'MarkerEdgeColor',[0 0.8 0]);
plot(PC1(celltypes{2}),u4(celltypes{2}),'o','MarkerSize',12,'MarkerFaceColor',[0.8 0 0],'MarkerEdgeColor',[0.8 0 0]);

keyboard

% Find putative ChATs
pinx = PC1 > 38;  % points clustered around ChATs
% pinx = find(PC1>27&PC1<36);
inxset = [ChATactinx PVactinx NTactinx];   % indices for conversion
rpinx = inxset(pinx);   % convert indices to ones relative to the PSTH matrix
pids = tags(rpinx);   % cellIDs
for k = 1:length(pids)   % open corresponding PSTH figures
    cid = pids{k};
    uiopen(['c:\Balazs\_analysis\NB\responsesorter3\' ...
        cid(1:end-2) '_' cid(end) '_PSTH.fig'],1)
end

% Overlay putative ChATs on LDA plot
plot(PC1(setdiff(find(pinx),celltypes{1})),u4(setdiff(find(pinx),celltypes{1})),...
    'o','MarkerSize',8,'MarkerFaceColor',[0 0.5 0],'MarkerEdgeColor',[0 0.5 0]);

% -------------------------------------------------------------------------
function C = tcorr(X,template)

% Correlation with specific PSTHs
C = template * X';