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
load([DATAPATH 'NB\responsesorter_new4\allPSTH.mat'])

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
selstr = ['"ChAT+"==0&"PV+"==0&"ID_PC">20&"Lr_PC"<0.15&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
NT = selectcell(selstr);   % cell IDs for ChAT cells

[nms inxa PVinx] = intersect(PV,tags);   % indices for PV cells
[nms inxa ChATinx] = intersect(ChAT,tags);   % indices for ChAT cells
[nms inxa NTinx] = intersect(NT,tags);   % indices for unidentified cells
PVactinx = intersect(PVinx,activated);   % indices for activated PV cells
ChATactinx = intersect(ChATinx,activated);   % indices for activated ChAT cells
NTactinx = intersect(NTinx,activated);  % indices for activated unidentified cells
numChAT = length(ChATactinx);   % number of ChAT cells
numPV = length(PVactinx);   % number of PV cells
numNT = length(NTactinx);   % number of unidentified cells
ChAT_PSTH = allpsth(ChATactinx,:);   % PSTHs of ChAT cells
PV_PSTH = allpsth(PVactinx,:);    % PSTHs of PV cells
NT_PSTH = allpsth(NTactinx,:);    % PSTHs of NT cells

% Correlation with specific PSTHs
template = mean(ChAT_PSTH);   % mean PSTH of ChAT cells
ChATcorr = ccorr(allpsth,template);   % correlation with mean ChAT PSTH
template = mean(NT_PSTH);   % mean PSTH of ChAT cells
NTcorr = ccorr(allpsth,template);   % correlation with mean PV PSTH

% Variables
v01 = activation_peak;   % peak time of effect
v02 = log10(maxvalue./baseline);   % relative firing rate change on log scale
v03 = activation_end - activation_start;   % activation duration
v04 = baseline;   % base firing rate
v05 = ChATcorr;   % correlation with mean ChAT PSTH
v06 = NTcorr;   % correlation with mean PV PSTH

% Restrict to ChAT and PV groups
v1 = [v01(ChATactinx) v01(NTactinx)]';
v2 = [v02(ChATactinx) v02(NTactinx)]';
v3 = [v03(ChATactinx) v03(NTactinx)]';
v4 = [v04(ChATactinx) v04(NTactinx)]';
v5 = [v05(ChATactinx) v05(NTactinx)]';
v6 = [v06(ChATactinx) v06(NTactinx)]';

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
vars = [v1 v3 v4 v5 v6];
celltyps = [zeros(numChAT,1); ones(numNT,1)];
W = lda2(vars,celltyps);   % linear discriminant analysis

% Add non-identified cells
u1 = [v01(ChATactinx) v01(NTactinx) v01(PVactinx)]';
u2 = [v02(ChATactinx) v02(NTactinx) v02(PVactinx)]';
u3 = [v03(ChATactinx) v03(NTactinx) v03(PVactinx)]';
u4 = [v04(ChATactinx) v04(NTactinx) v04(PVactinx)]';
u5 = [v05(ChATactinx) v05(NTactinx) v05(PVactinx)]';
u6 = [v06(ChATactinx) v06(NTactinx) v06(PVactinx)]';

% 'Standardize' new variables
u1 = (u1 - meanv1) / stdv1;   % standardize using the original group means and SDs
u2 = (u2 - meanv2) / stdv2;
u3 = (u3 - meanv3) / stdv3;
u4 = (u4 - meanv4) / stdv4;
u5 = (u5 - meanv5) / stdv5;
u6 = (u6 - meanv6) / stdv6;

% Plot LDA 
vars2 = [u1 u3 u4 u5 u6];
scores = [ones(size(vars2,1),1) vars2] * W';    % project all groups (including all non-tagged) to the subspace found
score_1 = scores(:,1);
score_2 = scores(:,2);
celltypes = {1:numChAT; numChAT+1:numChAT+numNT; numChAT+numNT+1:size(vars2,1)};

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
plot(PC1(celltypes{3}),PC2(celltypes{3}),'o','MarkerSize',5,'Color',[0.8 0.8 0.8],'MarkerFaceColor',[0.8 0.8 0.8]);
hold on
plot(PC1(celltypes{1}),PC2(celltypes{1}),'o','MarkerSize',12,'MarkerFaceColor',[0 0.8 0],'MarkerEdgeColor',[0 0.8 0]);
plot(PC1(celltypes{2}),PC2(celltypes{2}),'o','MarkerSize',5,'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7]);

% Plot LDA vs baseline firing rate
% u2 = [v02(ChATactinx) v02(NTactinx) v02(PVactinx)]';  % for non-normalized firing rate
figure   % plot
plot(PC1(celltypes{3}),u2(celltypes{3}),'o','MarkerSize',3,'Color','k','MarkerFaceColor','k');
hold on
plot(PC1(celltypes{1}),u2(celltypes{1}),'o','MarkerSize',8,'MarkerFaceColor',[0 0.8 0],'MarkerEdgeColor',[0 0.8 0]);
plot(PC1(celltypes{2}),u2(celltypes{2}),'o','MarkerSize',3,'MarkerFaceColor','k','MarkerEdgeColor','k');

keyboard

% Find putative ChATs
% pinx = PC1 > 38;  % points clustered around ChATs
pinx = find(PC1>1&PC1<1.5);
inxset = [ChATactinx NTactinx PVactinx];   % indices for conversion
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

keyboard

l1 = 'peak time of effect';
l2 = 'relative firing rate change on log scale';
l3 = 'activation duration';
l4 = 'base firing rate';
l5 = 'correlation with mean ChAT PSTH';
l6 = 'correlation with mean PV PSTH';

xa = 2;
xb = 5;
va = eval(['v0' num2str(xa)]);
vb = eval(['v0' num2str(xb)]);
la = eval(['l' num2str(xa)]);
lb = eval(['l' num2str(xb)]);

% va = v02 ./ v03;
% va = v02 ./ v03 ./ v04;

figure;plot(va(NTactinx),vb(NTactinx),'k.')
hold on;plot(va(ChATactinx),vb(ChATactinx),'go','MarkerSize',8,'MarkerFaceColor','g')
xlabel(la)
ylabel(lb)

figure
plot3(PC1(celltypes{2}),u4(celltypes{2}),u6(celltypes{2}),'o','MarkerSize',5,'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7]);
hold on
plot3(PC1(celltypes{1}),u4(celltypes{1}),u6(celltypes{1}),'o','MarkerSize',12,'MarkerFaceColor',[0 0.8 0],'MarkerEdgeColor',[0 0.8 0]);

% -------------------------------------------------------------------------
function C = tcorr(X,template)

% Correlation with specific PSTHs
C = template * X';

% -------------------------------------------------------------------------
function C = ccorr(X,template)

% Correlation with specific PSTHs
NumCells = size(X,1);
C = nan(NumCells);
for k = 1:NumCells
    xc = xcorr(template,X(k,:),10);
    C(k) = max(xc);
end