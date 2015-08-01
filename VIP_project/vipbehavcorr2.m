function vipbehavcorr2
%VIPBEHAVCORR2   Behavioral response profiles of VIP neurons.
%   VIPBEHAVCORR2 sorts PETHs calculated by VIPRESPONSESORTER into groups
%   based on significant inhibition/activation around behavioral feedback
%   and VIP positivity. Significance is determined based on cue
%   tone-aligned PETHs. Average stimulus-aligned PETHs are plotted.
%
%   See also VIPRESPONSESORTER.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   16-May-2013

% First CellBase: air-puff
choosecb('VIP_gonogo')
global DATAPATH
load([DATAPATH 'VIP\responsesorter2_Hit_newwin4_stimaligned\allPSTH.mat'])  % load PSTHs
[VIPactinx1 NTactinx1 inhibited1] = groups(allstats,tags);   % grouping of cells based on feedback-alignment
load([DATAPATH 'VIP\responsesorter2_Hit_newwin4\allPSTH.mat'])  % load PSTHs
allpsth1 = allpsth;   % PSTH
allspsth1 = smoothPSTH(allpsth1);   % smoothed PSTH

% Second CellBase: foot-shock
choosecb('VIP_gonogo2')
load([DATAPATH 'VIP\responsesorter2_Hit_newwin4_stimaligned_shock\allPSTH.mat'])  % load PSTHs
[VIPactinx2 NTactinx2 inhibited2] = groups(allstats,tags);   % grouping of cells based on feedback-alignment
load([DATAPATH 'VIP\responsesorter2_Hit_newwin4_shock\allPSTH.mat'])  % load PSTHs
allpsth2 = allpsth;   % PSTH
allspsth2 = smoothPSTH(allpsth2);   % smoothed PSTH

% Plot
plotPSTH(VIPactinx1,NTactinx1,inhibited1,allpsth1,allspsth1,...
    VIPactinx2,NTactinx2,inhibited2,allpsth2,allspsth2,time)

% -------------------------------------------------------------------------
function [VIPactinx NTactinx inhibited] = groups(allstats,tags)

% VIP+ and VIP- cells
selstr = '"VIP+"==1&"validity"==1';
VIP = selectcell(selstr);   % cellIDs of VIP+ cells
selstr = '"VIP+"~=1&"validity"==1';
NT = selectcell(selstr);   % cellIDs of VIP- cells

% Response profile variables
activation_peak = [allstats.activation_peak];   % peak time of activation
activation_start = [allstats.activation_start];   % start time of activation
activation_end = [allstats.activation_end];   % end time of activation
maxvalue = [allstats.maxvalue];   % maximal FR
Wpa = [allstats.Wpa];   % Mann-Whitney test for significant activation
inhibition_peak = [allstats.inhibition_peak];   % peak time of inhibition
inhibition_start = [allstats.inhibition_start];   % start time of inhibition
inhibition_end = [allstats.inhibition_end];   % end time of inhibition
minvalue = [allstats.minvalue];   % minimal FR
Wpi = [allstats.Wpi];   % Mann-Whitney test for significant inhibition
baseline = [allstats.baseline];   % baseline FR

% Groups of activated and inhibited cells
inx_act = Wpa < 0.01;   % significant activation
inx_inh = Wpi < 0.01;   % significant inhibition
vldty = getvalue('validity')';
activated = find(inx_act&~inx_inh&vldty);    % indices for activated cells
inhibited = find(inx_inh&~inx_act&vldty);    % indices for inhibited cells
ai = find(inx_act&inx_inh&vldty);   % indices for activated and inhibited cells
inx = activation_peak(ai) < inhibition_peak(ai);
activated_inhibited = sort(ai(inx));
inhibited_activated = ai(~inx);
activated = [activated'; activated_inhibited'];   % categorize based on first significant effect
inhibited = [inhibited'; inhibited_activated'];

[~, ~, VIPinx] = intersect(VIP,tags);   % indices for VIP cells
[~, ~, NTinx] = intersect(NT,tags);   % indices for VIP- cells
VIPactinx = intersect(VIPinx,activated);   % indices for activated VIP cells
NTactinx = setdiff(activated,VIPactinx);  % indices for unidentified cells
numVIP = length(VIPactinx);   % number of VIP cells
numNT = length(NTactinx);   % number of unidentified cells

% Population PSTHs for groups
% figure;imagesc(allpsth(NTactinx,:));colormap hot
% set(gca,'CLim',[0 5])
% figure;imagesc(allpsth(inhibited,:));colormap hot
% set(gca,'CLim',[0 5])
% figure;imagesc(allpsth(VIPactinx,:));colormap hot
% set(gca,'CLim',[0 5])

% -------------------------------------------------------------------------
function plotPSTH(VIPactinx1,NTactinx1,inhibited1,allpsth1,allspsth1,...
    VIPactinx2,NTactinx2,inhibited2,allpsth2,allspsth2,time)

% Colors
brown = [0.32 0.19 0.19];
purple = [0.48 0.06 0.89];
grey1 = [0.7 0.7 0.7];
grey2 = [0.4 0.4 0.4];
green = [0.2 0.8 0.2];
red = [0.8 0 0];
pink = [1 0.6 0.78];

% Average PSTH
% NTpsth = [allpsth1(NTactinx1,:); allpsth2(NTactinx2,:)];  % pooled NT PSTH matrix
% VIPpsth = [allpsth1(VIPactinx1,:); allpsth2(VIPactinx2,:)];  % pooled VIP PSTH matrix
% INHpsth = [allpsth1(inhibited1,:); allpsth2(inhibited2,:)];  % pooled inhibited PSTH matrix
% 
% figure
% hold on;
% errorshade(time,mean(NTpsth),std(NTpsth)/sqrt(size(NTpsth,1)),...
%     'LineColor',grey1,'ShadeColor',grey1)
% errorshade(time,mean(VIPpsth),std(VIPpsth)/sqrt(size(VIPpsth,1)),...
%     'LineColor',green,'ShadeColor',green)
% errorshade(time,mean(INHpsth),std(INHpsth)/sqrt(size(INHpsth,1)),...
%     'LineColor',pink,'ShadeColor',pink)
% line([0 0],[-4 20],'Color','k','LineStyle',':')
% line([-1200 2200],[0 0],'Color','k')
% set(gca,'box','off',...
%     'FontSize',16,'TickDir','out','XLim',[-200 600],'Ylim',[-2 4])
% xlabel('Time')
% ylabel('Normalized firing rate')

% Average PSTH, smoothed
NTpsth = [allspsth1(NTactinx1,:); allspsth2(NTactinx2,:)];  % pooled NT PSTH matrix
VIPpsth = [allspsth1(VIPactinx1,:); allspsth2(VIPactinx2,:)];  % pooled VIP PSTH matrix
INHpsth = [allspsth1(inhibited1,:); allspsth2(inhibited2,:)];  % pooled inhibited PSTH matrix

figure
hold on;
errorshade(time,mean(NTpsth),std(NTpsth)/sqrt(size(NTpsth,1)),...
    'LineColor',grey1,'ShadeColor',grey1)
errorshade(time,mean(VIPpsth),std(VIPpsth)/sqrt(size(VIPpsth,1)),...
    'LineColor',green,'ShadeColor',green)
errorshade(time,mean(INHpsth),std(INHpsth)/sqrt(size(INHpsth,1)),...
    'LineColor',pink,'ShadeColor',pink)
line([0 0],[-4 20],'Color','k','LineStyle',':')
line([-1200 2200],[0 0],'Color','k')
set(gca,'box','off',...
    'FontSize',16,'TickDir','out','XLim',[-200 600],'Ylim',[-2 4])
xlabel('Time')
ylabel('Normalized firing rate')

% Average PSTH, smoothed and baseline-corrected
VIPbaseline = VIPpsth(:,1:500);
NTbaseline = NTpsth(:,1:500);
INHbaseline = INHpsth(:,1:500);
figure
hold on;
errorshade(time,mean(NTpsth)-mean(NTbaseline(:)),std(NTpsth)/sqrt(size(NTpsth,1)),...
    'LineColor',grey1,'ShadeColor',grey1)
errorshade(time,mean(VIPpsth)-mean(VIPbaseline(:)),std(VIPpsth)/sqrt(size(VIPpsth,1)),...
    'LineColor',green,'ShadeColor',green)
errorshade(time,mean(INHpsth)-mean(INHbaseline(:)),std(INHpsth)/sqrt(size(INHpsth,1)),...
    'LineColor',pink,'ShadeColor',pink)
line([0 0],[-4 20],'Color','k','LineStyle',':')
line([-1200 2200],[0 0],'Color','k')
set(gca,'box','off',...
    'FontSize',16,'TickDir','out','XLim',[-1000 2000],'Ylim',[-3 5])
xlabel('Time')
ylabel('Normalized firing rate')

keyboard

% -------------------------------------------------------------------------
function spsth = smoothPSTH(psth)

NumCells = size(psth,1);   % number of cells
spsth = nan(size(psth));
for k = 1:NumCells   % loop through all cells
    spsth(k,:) = smooth(psth(k,:),'linear',21);   % smooth by moving average
end