function nbclustering
%NBCLUSTERING   Clustering of behavioral responses.
%   NBCLUSTERING performs hierarchical clustering (using Ward's
%   amalgamation rule) on parameters of PETHs aligned to behavioral events.
%   Baseline firing rate, relative firing rate change and correlation with
%   average response of identified cells are used as variables, extracted
%   from the PETHs aligned to false alarm responses (see NBRESPONSESORTER).
%   Only significantly activated (p<0.01, Mann-Whitney U-test; see
%   ULTIMATE_PSTH and PSTH_STATS) cells are included. Relative firing rate
%   change is calculated as logarithm of maximal firing rate divided by
%   baseline firing rate. Correlation is calculated as maximum of PETH
%   cross-correlation within a 10 ms window. All variables are z-scored.
%
%   See also NBRESPONSESORTER, ULTIMATE_PSTH, PSTH_STATS and NBLDA.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   31-Aug-2012

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
minvalue = [allstats.minvalue];   % minimal firing rate
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
selstr = ['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
ChAT = selectcell(selstr);   % cell IDs for ChAT cells
selstr = ['"ChAT+"==0&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
NT = selectcell(selstr);   % cell IDs for ChAT cells

[nms inxa ChATinx] = intersect(ChAT,tags);   % indices for ChAT cells
[nms inxa NTinx] = intersect(NT,tags);   % indices for unidentified cells
ChATactinx = intersect(ChATinx,activated);   % indices for activated ChAT cells
NTactinx = intersect(NTinx,activated);  % indices for activated unidentified cells
numChAT = length(ChATactinx);   % number of ChAT cells
numNT = length(NTactinx);   % number of unidentified cells
ChAT_PSTH = allpsth(ChATactinx,:);   % PSTHs of ChAT cells
NT_PSTH = allpsth(NTactinx,:);    % PSTHs of NT cells

% Correlation with specific PSTHs (cross-correlation-peak within 10 ms)
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

% Restrict to activated ChAT and NT groups
v1 = [v01(ChATactinx) v01(NTactinx)]';
v2 = [v02(ChATactinx) v02(NTactinx)]';
v3 = [v03(ChATactinx) v03(NTactinx)]';
v4 = [v04(ChATactinx) v04(NTactinx)]';
v5 = [v05(ChATactinx) v05(NTactinx)]';
v6 = [v06(ChATactinx) v06(NTactinx)]';
ChATvinx = 1:length(ChATactinx);   % new indices for ChAT cells
NTvinx = length(ChATactinx)+1:length(v1);   % new indices for NT cells

% Standardize variables
v1 = standardize(v1);
v2 = standardize(v2);
v3 = standardize(v3);
v4 = standardize(v4);
v5 = standardize(v5);
v6 = standardize(v6);

% Clustering
dmtx = [v4 v2 v5];
dist = pdist(dmtx);
links = linkage(dist,'Ward');
c = cluster(links,2);
pChATvinx = setdiff(find(c==1),ChATvinx);   % putative cholinergic neurons, new indexing
inxset = [ChATactinx NTactinx];   % indices for index conversion
pChATactinx = inxset(pChATvinx);   % putative cholinergic neurons, old indexing
pChAT = tags(pChATactinx);   % cell IDs of pChAT cells
NTvinx2 = setdiff(NTvinx,pChATvinx);   % remove pChAT from NT

% Gap statistic
% Ccc = cophenet(links,dist);
% [khat Gap s_k C] = gap_statistics(dmtx,20);

% Labels for plotting
l1 = 'peak time of effect';
l2 = 'relative firing rate change on log scale';
l3 = 'activation duration';
l4 = 'base firing rate';
l5 = 'correlation with mean ChAT PSTH';
l6 = 'correlation with mean PV PSTH';

% Colors
grey = [0.7 0.7 0.7];
green = [0 0.8 0];
blue = [0 0 0.7];

% 2D scatter plot of features
xa = 2;  % select feature for x axis
xb = 5;  % select feature for y axis
va = eval(['v0' num2str(xa)]);   % define x axis feature
vb = eval(['v0' num2str(xb)]);   % define y axis feature
la = eval(['l' num2str(xa)]);   % x axis label
lb = eval(['l' num2str(xb)]);   % y axis label
figure
plot(va(NTactinx),vb(NTactinx),'k.')   % NT cells
hold on
plot(va(ChATactinx),vb(ChATactinx),'o','MarkerSize',8,'MarkerFaceColor',...
    green,'MarkerEdgeColor',green)   % ChAT cells
plot(va(pChATactinx),vb(pChATactinx),'o','MarkerSize',8,'MarkerFaceColor',...
    blue,'MarkerEdgeColor',blue);   % pChAT cells
xlabel(la)
ylabel(lb)

% 3D scatter plots
figure
plot3(v1(NTvinx2),v2(NTvinx2),v5(NTvinx2),'o','MarkerSize',5,...
    'MarkerFaceColor', grey,'MarkerEdgeColor',grey);   % NT cells
hold on
plot3(v1(pChATvinx),v2(pChATvinx),v5(pChATvinx),'o','MarkerSize',12,...
    'MarkerFaceColor',blue,'MarkerEdgeColor',blue);   % pChAT cells
plot3(v1(ChATvinx),v2(ChATvinx),v5(ChATvinx),'o','MarkerSize',12,...
    'MarkerFaceColor',green,'MarkerEdgeColor',green);   % ChAT cells
grid on

figure
plot3(v3(NTvinx2),v2(NTvinx2),v5(NTvinx2),'o','MarkerSize',5,...
    'MarkerFaceColor', grey,'MarkerEdgeColor',grey);   % NT cells
hold on
plot3(v3(pChATvinx),v2(pChATvinx),v5(pChATvinx),'o','MarkerSize',12,...
    'MarkerFaceColor',blue,'MarkerEdgeColor',blue);   % pChAT cells
plot3(v3(ChATvinx),v2(ChATvinx),v5(ChATvinx),'o','MarkerSize',12,...
    'MarkerFaceColor',green,'MarkerEdgeColor',green);   % ChAT cells
grid on

figure
plot3(v4(NTvinx2),v2(NTvinx2),v5(NTvinx2),'o','MarkerSize',5,...
    'MarkerFaceColor', grey,'MarkerEdgeColor',grey);   % NT cells
hold on
plot3(v4(pChATvinx),v2(pChATvinx),v5(pChATvinx),'o','MarkerSize',12,...
    'MarkerFaceColor',blue,'MarkerEdgeColor',blue);   % pChAT cells
plot3(v4(ChATvinx),v2(ChATvinx),v5(ChATvinx),'o','MarkerSize',12,...
    'MarkerFaceColor',green,'MarkerEdgeColor',green);   % ChAT cells
grid on

keyboard

% -------------------------------------------------------------------------
function C = ccorr(X,template)

% Correlation with specific PSTHs
NumCells = size(X,1);
C = nan(1,NumCells);
for k = 1:NumCells
    xc = xcorr(template,X(k,:),10);
    C(k) = max(xc);
end