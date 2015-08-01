%% colors

brown = [0.32 0.19 0.19];
purple = [0.48 0.06 0.89];
grey1 = [0.7 0.7 0.7];
grey2 = [0.4 0.4 0.4];
green = [0 0.8 0];
red = [0.8 0 0];

%% load 

global DATAPATH
load([DATAPATH 'VIP\responsesorter_FA_newwin\allPSTH.mat'])


%% groups

selstr = '"VIP+"==1&"validity"==1';
VIP = selectcell(selstr);   % cellIDs of VIP+ cells
selstr = '"VIP+"~=1&"validity"==1';
NT = selectcell(selstr);   % cellIDs of VIP- cells

activation_peak_FA = [allstats.activation_peak];   % peak time of activation
activation_start_FA = [allstats.activation_start];
activation_end_FA = [allstats.activation_end];
maxvalue_FA = [allstats.maxvalue];
Wpa = [allstats.Wpa];   % Mann-Whitney test for significant activation
inhibition_peak_FA = [allstats.inhibition_peak];   % peak time of inhibition
inhibition_start_FA = [allstats.inhibition_start];
inhibition_end_FA = [allstats.inhibition_end];
minvalue_FA = [allstats.minvalue];
Wpi = [allstats.Wpi];   % Mann-Whitney test for significant inhibition
baseline_FA = [allstats.baseline];

% Groups of activated and inhibited cells
inx_act = Wpa < 0.01;   % significant activation
inx_inh = Wpi < 0.01;   % significant inhibition
vldty = getvalue('validity')';
activated = find(inx_act&~inx_inh&vldty);    % indices for activated cells
inhibited = find(inx_inh&~inx_act&vldty);    % indices for inhibited cells
ai = find(inx_act&inx_inh&vldty);   % indices for activated and inhibited cells
inx = activation_peak_FA(ai) < inhibition_peak_FA(ai);
activated_inhibited = sort(ai(inx));
inhibited_activated = ai(~inx);
activated_FA = [activated'; activated_inhibited'];   % categorize based on first significant effect
inhibited = [inhibited'; inhibited_activated'];

[nms inxa VIPinx] = intersect(VIP,tags);   % indices for VIP cells
VIPactinx_FA = intersect(VIPinx,activated);   % indices for activated VIP cells
NTactinx_FA = setdiff(activated,VIPactinx_FA);  % indices for unidentified cells
numVIP = length(VIPactinx_FA);   % number of VIP cells
numNT = length(NTactinx_FA);   % number of unidentified cells
VIP_PSTH = allpsth(VIPactinx_FA,:);    % PSTHs of VIP cells
[~, ~, NTinx] = intersect(NT,tags);   % indices for VIP- cells

% Light effects
lightinhinx = [2 18 25 30 59 63 73 75 83 86 88];
lightinhinx2 = [2 18 30 59 63 75 83 86 88];  % inh-act not included
lightactinx = [33 68 80 84 89];
lightactinx2 = [33 68 84 89];  % act-inh not included

%% load 

global DATAPATH
load([DATAPATH 'VIP\responsesorter_Hit_newwin\allPSTH.mat'])


%% groups

selstr = '"VIP+"==1&"validity"==1';
VIP = selectcell(selstr);   % cellIDs of VIP+ cells
selstr = '"VIP+"~=1&"validity"==1';
NT = selectcell(selstr);   % cellIDs of VIP- cells

activation_peak_Hit = [allstats.activation_peak];   % peak time of activation
activation_start_Hit = [allstats.activation_start];
activation_end_Hit = [allstats.activation_end];
maxvalue_Hit = [allstats.maxvalue];
Wpa = [allstats.Wpa];   % Mann-Whitney test for significant activation
inhibition_peak_Hit = [allstats.inhibition_peak];   % peak time of inhibition
inhibition_start_Hit = [allstats.inhibition_start];
inhibition_end_Hit = [allstats.inhibition_end];
minvalue_Hit = [allstats.minvalue];
Wpi = [allstats.Wpi];   % Mann-Whitney test for significant inhibition
baseline_Hit = [allstats.baseline];

% Groups of activated and inhibited cells
inx_act = Wpa < 0.01;   % significant activation
inx_inh = Wpi < 0.01;   % significant inhibition
vldty = getvalue('validity')';
activated = find(inx_act&~inx_inh&vldty);    % indices for activated cells
inhibited = find(inx_inh&~inx_act&vldty);    % indices for inhibited cells
ai = find(inx_act&inx_inh&vldty);   % indices for activated and inhibited cells
inx = activation_peak_Hit(ai) < inhibition_peak_Hit(ai);
activated_inhibited = sort(ai(inx));
inhibited_activated = ai(~inx);
activated_Hit = [activated'; activated_inhibited'];   % categorize based on first significant effect
inhibited = [inhibited'; inhibited_activated'];

%%

figure
plot(maxvalue_FA(NTinx)./baseline_FA(NTinx),maxvalue_Hit(NTinx)./baseline_Hit(NTinx),'k.')
hold on
plot(maxvalue_FA(VIPinx)./baseline_FA(VIPinx),maxvalue_Hit(VIPinx)./baseline_Hit(VIPinx),'go')

%%

figure
plot(maxvalue_FA(NTinx)-baseline_FA(NTinx),maxvalue_Hit(NTinx)-baseline_Hit(NTinx),'k.')
hold on
plot(maxvalue_FA(VIPinx)-baseline_FA(VIPinx),maxvalue_Hit(VIPinx)-baseline_Hit(VIPinx),'go')

%%

figure
plot(maxvalue_FA(NTinx),maxvalue_Hit(NTinx),'k.')
hold on
plot(maxvalue_FA(VIPinx),maxvalue_Hit(VIPinx),'go','MarkerEdgeColor',green,'MarkerFaceColor',green)

%%

figure
plot(maxvalue_FA(NTinx),activation_end_Hit(NTinx),'k.')
hold on
plot(maxvalue_FA(VIPinx),activation_end_Hit(VIPinx),'go','MarkerEdgeColor',green,'MarkerFaceColor',green)

%%

% figure
% plot(maxvalue_FA(NTactinx),maxvalue_Hit(NTactinx),'k.')
% hold on
% plot(maxvalue_FA(VIPactinx),maxvalue_Hit(VIPactinx),'go','MarkerEdgeColor',green,'MarkerFaceColor',green)


%% pdf #3 - accepted version

edges = 0:2:122;
centers = (edges(1:end-1) + edges(2:end)) / 2;
dist_FA_NT = histc(maxvalue_FA(NTinx),edges);
dist_FA_NT = dist_FA_NT(1:end-1);
dist_FA_NT = dist_FA_NT / sum(dist_FA_NT);

dist_FA_VIP = histc(maxvalue_FA(VIPinx),edges);
dist_FA_VIP = dist_FA_VIP(1:end-1);
dist_FA_VIP = dist_FA_VIP / sum(dist_FA_VIP);

dist_Hit_NT = histc(maxvalue_Hit(NTinx),edges);
dist_Hit_NT = dist_Hit_NT(1:end-1);
dist_Hit_NT = dist_Hit_NT / sum(dist_Hit_NT);

dist_Hit_VIP = histc(maxvalue_Hit(VIPinx),edges);
dist_Hit_VIP = dist_Hit_VIP(1:end-1);
dist_Hit_VIP = dist_Hit_VIP / sum(dist_Hit_VIP);

figure
plot(centers,dist_FA_NT,'Color','k','LineWidth',3)
hold on
plot(centers,dist_FA_VIP,'Color',green,'LineWidth',3)

figure
plot(centers,dist_Hit_NT,'Color','k','LineWidth',3)
hold on
plot(centers,dist_Hit_VIP,'Color',green,'LineWidth',3)

%% cdf #3

edges = -2:2:122;
centers = (edges(1:end-1) + edges(2:end)) / 2;
dist_FA_NT = histc(maxvalue_FA(NTinx),edges);
dist_FA_NT = dist_FA_NT(1:end-1);
dist_FA_NT = dist_FA_NT / sum(dist_FA_NT);

dist_FA_VIP = histc(maxvalue_FA(VIPinx),edges);
dist_FA_VIP = dist_FA_VIP(1:end-1);
dist_FA_VIP = dist_FA_VIP / sum(dist_FA_VIP);

dist_Hit_NT = histc(maxvalue_Hit(NTinx),edges);
dist_Hit_NT = dist_Hit_NT(1:end-1);
dist_Hit_NT = dist_Hit_NT / sum(dist_Hit_NT);

dist_Hit_VIP = histc(maxvalue_Hit(VIPinx),edges);
dist_Hit_VIP = dist_Hit_VIP(1:end-1);
dist_Hit_VIP = dist_Hit_VIP / sum(dist_Hit_VIP);

figure
plot(centers,cumsum(dist_FA_NT),'k','LineWidth',3)
hold on
plot(centers,cumsum(dist_FA_VIP),'Color',green,'LineWidth',3)
xlim([edges(1) edges(end)])
ylim([0 1.02])
box off

figure
plot(centers,cumsum(dist_Hit_NT),'k','LineWidth',3)
hold on
plot(centers,cumsum(dist_Hit_VIP),'Color',green,'LineWidth',3)
xlim([edges(1) edges(end)])
ylim([0 1.02])
box off