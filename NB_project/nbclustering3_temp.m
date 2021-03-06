function nbclustering3
%NBCLUSTERING2   Clustering of behavioral responses.
%   NBCLUSTERING2 performs hierarchical clustering (using Ward's
%   amalgamation rule) on parameters of PETHs aligned to behavioral events.
%   Baseline firing rate, relative firing rate change and correlation with
%   average response of identified cells are used as variables, extracted
%   from the PETHs aligned to false alarm responses (see NBRESPONSESORTER).
%   Only significantly activated (p<0.01, Mann-Whitney U-test; see
%   ULTIMATE_PSTH and PSTH_STATS) cells are included. Relative firing rate
%   change is calculated as logarithm of maximal firing rate divided by
%   baseline firing rate. Correlation is calculated as maximum of PETH
%   cross-correlation within an 8 ms window. All variables are z-scored.
%
%   See also NBRESPONSESORTER, ULTIMATE_PSTH, PSTH_STATS and NBLDA.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   31-Aug-2012

% Edit log: BH 8/31/12

% Directories
global DATAPATH
resdir = ([DATAPATH 'NB\clustering_newdata\']);   % result directory

% Cells
selstr = ['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
ChAT = selectcell(selstr);   % cell IDs for ChAT cells
ChAT = [ChAT 'n045_121217x_4.6'];   % clustered based on light-evoked spikes
selstr = ['"pChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
pChAT = selectcell(selstr);  % putative
selstr = ['"ChAT+"==0&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
NT = selectcell(selstr);   % cell IDs for non-tagged cells (note: some 
% pChATs are included with more than one session here; they can only be
% checked for repetitions AFTER we identify them as pChATs)

% Load data
% tags = [ChAT NT];
% allpsth = getvalue('FA_psth',tags);
% allpsth = nancell2mat(allpsth);
% allpsth = zscore(allpsth,0,2);
% allstats = getvalue('FA_psth_stats',tags);
% allstats = nancell2struct2(allstats);

load([DATAPATH 'NB\responsesorter_new4\allPSTH.mat'])
ta = load([DATAPATH 'NB\responsesorter_n045_121217x_4_6\allPSTH.mat']);  % add the cell cluster based on light-evoked spike shape
ta2 = load([DATAPATH 'NB\responsesorter_newdata2\allPSTH.mat']);  % add the new cells
allstats = [allstats; ta.allstats; ta2.allstats]; %#ok<*NODEF>
allpsth = [allpsth; ta.allpsth; ta2.allpsth];
tags = [tags ta.tags ta2.tags];

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
% NOTE: 405 (not 409) activated cells because some pChATs were recorded
% more than once!

% Activated NB cells
[nms inxa ChATinx] = intersect(ChAT,tags);   % indices for ChAT cells
[nms inxa NTinx] = intersect(NT,tags);   % indices for unidentified cells
ChATactinx = intersect(ChATinx,activated);   % indices for activated ChAT cells
NTactinx = intersect(NTinx,[activated; inhibited]);  % indices for activated unidentified cells
numChAT = length(ChATactinx);   % number of ChAT cells
numNT = length(NTactinx);   % number of unidentified cells
ChAT_PSTH = allpsth(ChATactinx,:);   % PSTHs of ChAT cells
NT_PSTH = allpsth(NTactinx,:);    % PSTHs of NT cells
act_PSTH = allpsth(activated,:);    % PSTHs of activated cells
inh_PSTH = allpsth(inhibited,:);    % PSTHs of inhibited cells

% Correlation with specific PSTHs (cross-correlation-peak within 10 ms)
template = mean(ChAT_PSTH);   % mean PSTH of ChAT cells
ChATcorr = ccorr(allpsth,template);   % correlation with mean ChAT PSTH
template = mean(NT_PSTH);   % mean PSTH of ChAT cells
NTcorr = ccorr(allpsth,template);   % correlation with mean PV PSTH
template = mean(act_PSTH);   % mean PSTH of activated cells
ACTcorr = ccorr(allpsth,template);   % correlation with mean activated PSTH
template = mean(inh_PSTH);   % mean PSTH of inhibited cells
INHcorr = ccorr(allpsth,template);   % correlation with mean inhibited PSTH

% Variables
v01 = activation_peak;   % peak time of effect
v02 = log10(maxvalue./baseline);   % relative firing rate change on log scale
v03 = activation_end - activation_start;   % activation duration
v04 = baseline;   % base firing rate
v05 = ChATcorr;   % correlation with mean ChAT PSTH
v06 = NTcorr;   % correlation with mean PV PSTH
v07 = ACTcorr;
v08 = INHcorr;

% Restrict to activated ChAT and NT groups
v1 = [v01(ChATactinx) v01(NTactinx)]';
v2 = [v02(ChATactinx) v02(NTactinx)]';
v3 = [v03(ChATactinx) v03(NTactinx)]';
v4 = [v04(ChATactinx) v04(NTactinx)]';
v5 = [v05(ChATactinx) v05(NTactinx)]';
v6 = [v06(ChATactinx) v06(NTactinx)]';
v7 = [v07(ChATactinx) v07(NTactinx)]';
v8 = [v08(ChATactinx) v08(NTactinx)]';
ChATvinx = 1:length(ChATactinx);   % new indices for ChAT cells
NTvinx = length(ChATactinx)+1:length(v1);   % new indices for NT cells

% Standardize variables
v1 = standardize(nan2zero(v1));
v2 = standardize(v2);
v3 = standardize(nan2zero(v3));
v4 = standardize(v4);
v5 = standardize(v5);
v6 = standardize(v6);
v7 = standardize(v7);
v8 = standardize(v8);

% Clustering
% dmtx = [v7 v8 v5 v6 v2];
% dmtx = [v7 v8 v5 v2];
% dmtx = [v7 v8 v5 v6 v2 v3];
% dmtx = [v7 v8 v5 v3 v2];
dmtx = [v2 v3 v5 v8];   % not bad
dmtx = [v2 v4 v5];
dmtx = [v2 v3 v5];
dmtx = [v2 v5 v6 v7];   % good
dmtx = [v2 v5 v6 v8];   % best
dist = pdist(dmtx);
links = linkage(dist,'Ward');
c = cluster(links,4);
pChATvinx = setdiff(find(c==1),ChATvinx);   % putative cholinergic neurons, new indexing
inxset = [ChATactinx; NTactinx];   % indices for index conversion
pChATactinx = inxset(pChATvinx);   % putative cholinergic neurons, old indexing
pChAT = tags(pChATactinx);   % cell IDs of pChAT cells
NTvinx2 = setdiff(NTvinx,pChATvinx);   % remove pChAT from NT
C1inx = sort(inxset(c==1));
C2inx = sort(inxset(c==2));
C3inx = sort(inxset(c==3));
C4inx = sort(inxset(c==4));

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
l7 = 'correlation with mean ACT PSTH';
l8 = 'correlation with mean INH PSTH';

% Colors
grey = [0.7 0.7 0.7];
green = [0 0.8 0];
blue = [0 0 0.7];

% 2D scatter plot of features
xa = 5;  % select feature for x axis
xb = 8;  % select feature for y axis
va = eval(['v0' num2str(xa)]);   % define x axis feature
vb = eval(['v0' num2str(xb)]);   % define y axis feature
la = eval(['l' num2str(xa)]);   % x axis label
lb = eval(['l' num2str(xb)]);   % y axis label
figure
plot(va(NTactinx),vb(NTactinx),'k.')   % NT cells
hold on
plot(va(C4inx),vb(C4inx),'o','MarkerSize',6,'MarkerFaceColor',...
    [0 0.2 0.4],'MarkerEdgeColor',[0 0.2 0.4]);
plot(va(C1inx),vb(C1inx),'o','MarkerSize',6,'MarkerFaceColor',...
    blue,'MarkerEdgeColor',blue);
plot(va(C2inx),vb(C2inx),'o','MarkerSize',6,'MarkerFaceColor',...
    [0.7 0.6 0.3],'MarkerEdgeColor',[0.7 0.6 0.3]);
plot(va(C3inx),vb(C3inx),'o','MarkerSize',6,'MarkerFaceColor',...
    [0.7 0.2 0.4],'MarkerEdgeColor',[0.7 0.2 0.4]);
plot(va(ChATactinx),vb(ChATactinx),'o','MarkerSize',8,'MarkerFaceColor',...
    green,'MarkerEdgeColor',green)   % ChAT cells
% plot(va(pChATactinx),vb(pChATactinx),'o','MarkerSize',8,'MarkerFaceColor',...
%     blue,'MarkerEdgeColor',blue);   % pChAT cells
set(gca,'box','off','LineWidth',1)
xlabel(la)
ylabel(lb)
% saveas(gcf,fullfile(resdir,'ChATcorr_vs_FRchange.fig'))

% 3D scatter plots
figure
plot3(v7(NTvinx2),v8(NTvinx2),v5(NTvinx2),'o','MarkerSize',5,...
    'MarkerFaceColor', grey,'MarkerEdgeColor',grey);   % NT cells
hold on
plot3(v7(pChATvinx),v8(pChATvinx),v5(pChATvinx),'o','MarkerSize',12,...
    'MarkerFaceColor',blue,'MarkerEdgeColor',blue);   % pChAT cells
plot3(v7(ChATvinx),v8(ChATvinx),v5(ChATvinx),'o','MarkerSize',12,...
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
% saveas(gcf,fullfile(resdir,'threeD_scatter.fig'))

keyboard

% Bar graphs
x1 = [v01(ChATactinx) v01(NTactinx)];   % peak time
x = [x1(ChATvinx) x1(pChATvinx)];
y = x1(NTvinx2);
boxstat(x,y,'cholinergic','other')

x1 = [v02(ChATactinx) v02(NTactinx)];   % response magnitude
x = [x1(ChATvinx) x1(pChATvinx)];
y = x1(NTvinx2);
boxstat(x,y,'cholinergic','other')

x1 = [v03(ChATactinx) v03(NTactinx)];   % response duration
x = [x1(ChATvinx) x1(pChATvinx)];
y = x1(NTvinx2);
boxstat(x,y,'cholinergic','other')

x1 = [v04(ChATactinx) v04(NTactinx)];   % baseline FR
x = [x1(ChATvinx) x1(pChATvinx)];
y = x1(NTvinx2);
boxstat(x,y,'cholinergic','other')

% -------------------------------------------------------------------------
function C = ccorr(X,template)

% Correlation with specific PSTHs
NumCells = size(X,1);
C = nan(1,NumCells);
for k = 1:NumCells
    xc = xcorr(template,X(k,:),8);
    C(k) = max(xc);
end

% -------------------------------------------------------------------------
function [H Wp] = boxstat(d1,d2,l1,l2,alpha,str)
%BOXSTAT   Box plot with statistics.
%   [H WP] = BOXSTAT(D1,D2,L1,L2,ALPHA) plots box plot (H, figure handle)
%   of data sets D1 and D2 using x axis labels L1 and L2. It performs
%   Mann-Whitney U-test with significance level ALPHA (default is 0.05) and
%   returns p value (WP).
%
%   [H WP] = BOXSTAT(D1,D2,L1,L2,ALPHA,'PAIRED') performes a paired
%   non-parametric test (Wilcoxon's signed rank test) instead of the
%   Mann-Whitney U-test.
%
%   See also BOXPLOT.

% Input argument check
error(nargchk(4,6,nargin))
if nargin < 6
    str = 'nonpaired';
end
if nargin < 5 || isempty(alpha)
    alpha = 0.05;   % default significance level
end
d1 = d1(:)';   % row vector format
d2 = d2(:)';

% Box plot
H = figure;
boxplot([d1 d2],[zeros(size(d1)) ones(size(d2))],'labels',[{l1} {l2}],...
    'widths',[0.7 0.7]);
switch str
    case 'nonpaired'
        [Wp, Wh] = ranksum(d1,d2,'alpha',alpha);   % Mann-Whitney U-test
    case 'paired'
        [Wp, Wh] = signrank(d1,d2,'alpha',alpha);   % Wilcoxon signed rank test
    otherwise
        error('boxstat:inputArg','Unsupported input argument.')
end
if Wh
    clr = 'red';
else
    clr = 'black';
end
y_lim = ylim;
x_lim = xlim;
tpos1 = x_lim(1) + (x_lim(2) - x_lim(1)) / 5;
tpos2 = y_lim(1) + (y_lim(2)-y_lim(1)) * 4 / 5;
text(tpos1,tpos2,num2str(Wp),'Color',clr,'Horizontalalignment','center')