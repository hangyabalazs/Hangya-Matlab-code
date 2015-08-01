function nblda
%SOMSUMMARY_FCN   Linear discriminant analysis.
%   SOMSUMMARY_FCN calculates LDA (linear discriminant analysis) on PSTHs
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

% Load data
global DATAPATH
load([DATAPATH 'NB\responsesorter3\allPSTH.mat'])

% PSTH statistics
activation_peak = [allstats.activation_peak];
activation_start = [allstats.activation_start];
activation_end = [allstats.activation_end];
maxvalue = [allstats.maxvalue];
Wpa = [allstats.Wpa];
inhibition_peak = [allstats.inhibition_peak];
inhibition_start = [allstats.inhibition_start];
inhibition_end = [allstats.inhibition_end];
minvalue = [allstats.minvalue];
Wpi = [allstats.Wpi];
baseline = [allstats.baseline];

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
NTactinx = setdiff(activated,union(ChATactinx,PVactinx));
numChAT = length(ChATactinx);
numPV = length(PVactinx);
numNT = length(NTactinx);

% Variables
v01 = activation_peak;   % peak time of effect
v02 = log10(maxvalue./baseline);   % relative firing rate change on log scale
v03 = activation_end - activation_start;   % activation duration
v04 = baseline;   % base firing rate

% Restrict to ChAT and PV groups
v1 = [v01(ChATactinx) v01(PVactinx)]';
v2 = [v02(ChATactinx) v02(PVactinx)]';
v3 = [v03(ChATactinx) v03(PVactinx)]';
v4 = [v04(ChATactinx) v04(PVactinx)]';

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

% LDA
vars = [v1 v2 v3 v4];
celltyps = [zeros(numChAT,1); ones(numPV,1)];
W = lda2(vars,celltyps);   % linear discriminant analysis

% Add non-identified cells
u1 = [v01(ChATactinx) v01(PVactinx) v01(NTactinx)]';
u2 = [v02(ChATactinx) v02(PVactinx) v02(NTactinx)]';
u3 = [v03(ChATactinx) v03(PVactinx) v03(NTactinx)]';
u4 = [v04(ChATactinx) v04(PVactinx) v04(NTactinx)]';

% 'Standardize' new variables
u1 = (u1 - meanv1) / stdv1;   % standardize using the original group means and SDs
u2 = (u2 - meanv2) / stdv2;
u3 = (u3 - meanv3) / stdv3;
u4 = (u4 - meanv4) / stdv4;

% Plot LDA 
vars2 = [u1 u2 u3 u4];
scores = [ones(size(vars2,1),1) vars2] * W';    % project all groups (including all non-tagged) to the plane found
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
PC1 = PCAscores(:,1);
PC2 = PCAscores(:,2);

figure   % plot
plot(PC1(celltypes{3}),PC2(celltypes{3}),'o','MarkerSize',5,'Color',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7]);
hold on
plot(PC1(celltypes{1}),PC2(celltypes{1}),'o','MarkerSize',12,'MarkerFaceColor',[0 0.8 0],'MarkerEdgeColor',[0 0.8 0]);
plot(PC1(celltypes{2}),PC2(celltypes{2}),'o','MarkerSize',12,'MarkerFaceColor',[0.8 0 0],'MarkerEdgeColor',[0.8 0 0]);