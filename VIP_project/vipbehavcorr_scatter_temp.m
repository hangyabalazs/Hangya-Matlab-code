function vipbehavcorr_scatter
%VIPBEHAVCORR_SCATTER   Behavioral response of VIP neurons.
%   VIPBEHAVCORR_SCATTER plots maximal firing rate for false alarm vs. hit
%   trials based on PSTHs calculated by VIPRESPONSESORTER.
%
%   See also VIPBEHAVCORR and VIPRESPONSESORTER.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   17-May-2013

% First CellBase: air-puff
choosecb('VIP_gonogo')
global DATAPATH
load([DATAPATH 'VIP\responsesorter2_Hit_newwin4\allPSTH.mat'])  % load Hit PSTHs
[~, ~ , ~, maxvalue1_Hit] = groups(allstats,tags);   % grouping of cells based on response

load([DATAPATH 'VIP\responsesorter2_FA_newwin4\allPSTH.mat'])  % load FA PSTHs
[VIPactinx1 NTactinx1 inhibited1 maxvalue1_FA] = groups(allstats,tags);   % grouping of cells based on response

% Second CellBase: foot-shock
choosecb('VIP_gonogo2')
load([DATAPATH 'VIP\responsesorter2_Hit_newwin4_shock\allPSTH.mat'])  % load PSTHs
[~, ~, ~, maxvalue2_Hit] = groups(allstats,tags);   % grouping of cells based on response

load([DATAPATH 'VIP\responsesorter2_FA_newwin4_shock\allPSTH.mat'])  % load PSTHs
[VIPactinx2 NTactinx2 inhibited2 maxvalue2_FA] = groups(allstats,tags);   % grouping of cells based on response

% Load firing rates between tone and feedback, and after feedback
load([DATAPATH 'VIP\FRs\tone_and_fb_rates.mat'])

% Maximal firing rates
mx_NT_FA = [fr_FA_tone1(NTactinx1) fr_FA_tone2(NTactinx2)];  % max FR after FA, NT cells
mx_NT_Hit = [fr_Hit_tone1(NTactinx1) fr_Hit_tone2(NTactinx2)];  % max FR after Hit, NT cells
mx_VIP_FA = [fr_FA_tone1(VIPactinx1) fr_FA_tone2(VIPactinx2)];  % max FR after FA, VIP cells
mx_VIP_Hit = [fr_Hit_tone1(VIPactinx1) fr_Hit_tone2(VIPactinx2)];  % max FR after Hit, VIP cells

% Scatter plot
green = [0.2 0.8 0.2];
figure
plot(mx_NT_FA,mx_NT_Hit,'k.')
hold on
plot(mx_VIP_FA,mx_VIP_Hit,'go','MarkerEdgeColor',green,'MarkerFaceColor',green)
xlabel('False Alarm')
ylabel('Hit')

% Marginal CDFs
edges = -2:2:52;  % histogram bin edges
centers = (edges(1:end-1) + edges(2:end)) / 2;   % histogram bin centers
dist_FA_NT = histc(mx_NT_FA,edges);   % distribution of max FR values after FAs, NT cells
dist_FA_NT = dist_FA_NT(1:end-1);
dist_FA_NT = dist_FA_NT / sum(dist_FA_NT);   % normalize

dist_FA_VIP = histc(mx_VIP_FA,edges);   % distribution of max FR values after FAs, VIP cells
dist_FA_VIP = dist_FA_VIP(1:end-1);
dist_FA_VIP = dist_FA_VIP / sum(dist_FA_VIP);   % normalize

dist_Hit_NT = histc(mx_NT_Hit,edges);   % distribution of max FR values after Hits, NT cells
dist_Hit_NT = dist_Hit_NT(1:end-1);
dist_Hit_NT = dist_Hit_NT / sum(dist_Hit_NT);   % normalize

dist_Hit_VIP = histc(mx_VIP_Hit,edges);   % distribution of max FR values after Hits, VIP cells
dist_Hit_VIP = dist_Hit_VIP(1:end-1);
dist_Hit_VIP = dist_Hit_VIP / sum(dist_Hit_VIP);   % normalize

figure
stairs(centers,cumsum(dist_FA_NT),'k','LineWidth',3)
hold on
stairs(centers,cumsum(dist_FA_VIP),'Color',green,'LineWidth',3)
xlim([edges(1) edges(end)])
ylim([0 1.02])
box off
title('False Alarm')

figure
stairs(centers,cumsum(dist_Hit_NT),'k','LineWidth',3)
hold on
stairs(centers,cumsum(dist_Hit_VIP),'Color',green,'LineWidth',3)
xlim([edges(1) edges(end)])
ylim([0 1.02])
box off
title('Hit')

keyboard

% -------------------------------------------------------------------------
function [VIPactinx NTactinx inhibited maxvalue] = groups(allstats,tags)

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