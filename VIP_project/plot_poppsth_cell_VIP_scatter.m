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
VIPactinx = intersect(VIPinx,activated);   % indices for activated VIP cells
NTactinx = setdiff(activated,VIPactinx);  % indices for unidentified cells
numVIP = length(VIPactinx);   % number of VIP cells
numNT = length(NTactinx);   % number of unidentified cells
VIP_PSTH = allpsth(VIPactinx,:);    % PSTHs of VIP cells
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

figure
plot(maxvalue_FA(NTactinx),maxvalue_Hit(NTactinx),'k.')
hold on
plot(maxvalue_FA(VIPactinx),maxvalue_Hit(VIPactinx),'go','MarkerEdgeColor',green,'MarkerFaceColor',green)


%%

%% pdf #3 - accepted version

edges = -1:2:121;
centers = (edges(1:end-1) + edges(2:end)) / 2;
dist_activated = histc(activation_peak(activated),edges);
dist_activated = dist_activated(1:end-1);
dist_activated = dist_activated / sum(dist_activated);

dist_inhibited = histc(inhibition_peak(inhibited),edges);
dist_inhibited = dist_inhibited(1:end-1);
dist_inhibited = dist_inhibited / sum(dist_inhibited);

dist_tagged = histc(activation_peak(tagged),edges);
dist_tagged = dist_tagged(1:end-1);
dist_tagged = dist_tagged / sum(dist_tagged);

figure
plot(centers,dist_activated,'Color','r','LineWidth',3)
hold on
plot(centers,dist_inhibited,'Color','b','LineWidth',3)
plot(centers,dist_tagged,'Color','k','LineWidth',3)
plot(centers,dist_tagged,'Color','k','LineWidth',3)

%% cdf #3

% figure
% plot(centers,cumsum(dist_activated),'r','LineWidth',3)
% hold on
% plot(centers,cumsum(dist_inhibited),'b','LineWidth',3)
% plot(centers,cumsum(dist_tagged),'k','LineWidth',3)

% figure
% plot(sort(activation_peak(activated),'ascend'),cumsum(ones(1,length(activated)))/length(activated),'r','LineWidth',3)
% hold on
% plot(sort(inhibition_peak(inhibited),'ascend'),cumsum(ones(1,length(inhibited)))/length(inhibited),'b','LineWidth',3)
% plot(sort(activation_peak(tagged),'ascend'),cumsum(ones(1,length(tagged)))/length(tagged),'k','LineWidth',3)

edges = -0.5:1:100.5;
centers = (edges(1:end-1) + edges(2:end)) / 2;
dist_activated = histc(activation_peak(activated),edges);
dist_activated = dist_activated(1:end-1);
dist_activated = dist_activated / sum(dist_activated);

dist_inhibited = histc(inhibition_peak(inhibited),edges);
dist_inhibited = dist_inhibited(1:end-1);
dist_inhibited = dist_inhibited / sum(dist_inhibited);

dist_tagged = histc(activation_peak(tagged),edges);
dist_tagged = dist_tagged(1:end-1);
dist_tagged = dist_tagged / sum(dist_tagged);

figure
plot(centers,cumsum(dist_activated),'r','LineWidth',3)
hold on
plot(centers,cumsum(dist_inhibited),'b','LineWidth',3)
plot(centers,cumsum(dist_tagged),'k','LineWidth',3)