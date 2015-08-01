%% NB

% Identified NB cells
choosecb('NB')
selstr = ['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
ChAT = selectcell(selstr);   % cell IDs for ChAT cells
ChAT = [ChAT 'n045_121217x_4.6'];   % clustered based on light-evoked spikes
selstr = ['"ChAT+"==0&"pChAT+"==0&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
NT = selectcell(selstr);   % cell IDs for non-tagged cells (includes PV)

%% spike properties

SpikeShape_NT1 = getvalue('SpikeShape',NT);
SpikeShape_NT1 = nancell2struct(SpikeShape_NT1);
SpikeWidth_NT1 = [SpikeShape_NT1.SpikeWidth];
Spike_NT1 = {SpikeShape_NT1.Spike};
Spike_NT1 = cell2mat(cellfun(@(s)zscore(s'),Spike_NT1,'UniformOutput',false))';

SpikeShape_ChAT1 = getvalue('SpikeShape',ChAT);
SpikeShape_ChAT1 = nancell2struct(SpikeShape_ChAT1);
SpikeWidth_ChAT1 = [SpikeShape_ChAT1.SpikeWidth];
Spike_ChAT1 = {SpikeShape_ChAT1.Spike};
Spike_ChAT1 = cell2mat(cellfun(@(s)zscore(s'),Spike_ChAT1,'UniformOutput',false))';

%% HDB

% Identified HDB cells
choosecb('HDB')
selstr = ['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
    'ismember("Area1",{''HDB'',''SI'',''VP''})'];
ChAT = selectcell(selstr);   % cell IDs for ChAT cells
ChAT = [ChAT 'n067_141017a_1.3'];   % without miss-triggering it passes the cluster criteria
ChAT = [ChAT 'n067_141019a_5.2'];   % light spike assisted clustering
selstr = ['"ChAT+"==0&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
    'ismember("Area1",{''HDB'',''SI'',''VP''})'];
NT = selectcell(selstr);   % cell IDs for non-tagged cells

%% spike properties

SpikeShape_NT2 = getvalue('SpikeShape',NT);
SpikeShape_NT2 = nancell2struct(SpikeShape_NT2);
SpikeWidth_NT2 = [SpikeShape_NT2.SpikeWidth];
Spike_NT2 = {SpikeShape_NT2.Spike};
Spike_NT2 = cell2mat(cellfun(@(s)zscore(s'),Spike_NT2,'UniformOutput',false))';

SpikeShape_ChAT2 = getvalue('SpikeShape',ChAT);
SpikeShape_ChAT2 = nancell2struct(SpikeShape_ChAT2);
SpikeWidth_ChAT2 = [SpikeShape_ChAT2.SpikeWidth];
Spike_ChAT2 = {SpikeShape_ChAT2.Spike};
Spike_ChAT2 = cell2mat(cellfun(@(s)zscore(s'),Spike_ChAT2,'UniformOutput',false))';

%% combined spike properties

SpikeShape_NT = [SpikeShape_NT1 SpikeShape_NT2];
SpikeWidth_NT = [SpikeWidth_NT1 SpikeWidth_NT2];
Spike_NT = [Spike_NT1; Spike_NT2];

SpikeShape_ChAT = [SpikeShape_ChAT1 SpikeShape_ChAT2];
SpikeWidth_ChAT = [SpikeWidth_ChAT1 SpikeWidth_ChAT2];
Spike_ChAT = [Spike_ChAT1; Spike_ChAT2];

%% colors

green = [0 0.8 0];
blue = [0 0 0.8];
grey = [0.7 0.7 0.7];

%%

% Scatter plot
feature1 = 'SpikeWidth';
feature2 = 'PeakToPostValleyTime';
figure
plot([SpikeShape_NT.(feature1)],[SpikeShape_NT.(feature2)],'o','MarkerSize',2,...
    'MarkerEdgeColor',grey,'MarkerFaceColor',grey)
hold on
plot([SpikeShape_ChAT.(feature1)],[SpikeShape_ChAT.(feature2)],'o',...
    'MarkerEdgeColor',green,'MarkerFaceColor',green)
plot([SpikeShape_pChAT.(feature1)],[SpikeShape_pChAT.(feature2)],'o',...
    'MarkerEdgeColor',blue,'MarkerFaceColor',blue)
xlabel(feature1)
ylabel(feature2)

% Stats
median([SpikeShape_NT.(feature1)])
median([SpikeShape_ChAT.(feature1)])
boxstat([SpikeShape_NT.(feature1)],[SpikeShape_ChAT.(feature1)],'NT','ChAT')
title(feature1)

median([SpikeShape_NT.(feature2)])
median([SpikeShape_ChAT.(feature2)])
boxstat([SpikeShape_NT.(feature2)],[SpikeShape_ChAT.(feature2)],'NT','ChAT')
title(feature2)

figure
hold on
bar([1 2],[median([SpikeShape_NT.(feature2)])...
    median([SpikeShape_ChAT.(feature2)])],'EdgeColor','k','FaceColor','none')
errorbar([1 2],[median([SpikeShape_NT.(feature2)])...
    median([SpikeShape_ChAT.(feature2)])],...
    [se_of_median([SpikeShape_NT.(feature2)])...
    se_of_median([SpikeShape_ChAT.(feature2)])],'k+')
title(feature2)

%% more stats

median([SpikeShape_ChAT.(feature2)])
median([SpikeShape_pChAT.(feature2)])
boxstat([SpikeShape_ChAT.(feature2)],[SpikeShape_pChAT.(feature2)],'ChAT','pChAT')
title(feature2)

median([SpikeShape_NT.(feature2)])
median([SpikeShape_pChAT.(feature2)])
boxstat([SpikeShape_NT.(feature2)],[SpikeShape_pChAT.(feature2)],'NT','pChAT')
title(feature2)

%% average spike shape

figure
hold on
errorshade(1:32,mean(Spike_ChAT),nanse(Spike_ChAT),'LineColor',green,...
    'ShadeColor',green)
errorshade(1:32,mean(Spike_NT),nanse(Spike_NT),'LineColor',[0 0 0],...
    'ShadeColor',[0 0 0])


%% average spike shape w controling for depth

figure
hold on
errorshade(1:32,mean(Spike_ChAT(depth_ChAT>4200,:)),...
    nanse(Spike_ChAT(depth_ChAT>4200,:)),'LineColor',green,...
    'ShadeColor',green)
errorshade(1:32,mean(Spike_pChAT(depth_pChAT>4200,:)),...
    nanse(Spike_pChAT(depth_pChAT>4200,:)),'LineColor',blue,...
    'ShadeColor',blue)
errorshade(1:32,mean(Spike_NT),nanse(Spike_NT),'LineColor',[0 0 0],...
    'ShadeColor',[0 0 0])
errorshade(1:32,mean(Spike_PV),nanse(Spike_PV),'LineColor',[0.8 0 0],...
    'ShadeColor',[0.8 0 0])

%% tonic/phasic

% % Load auto-correlations
% load('C:\Balazs\_analysis\NB\ACG\medium\ACG_matrices.mat')

% ACG properties
Refractory_ChAT = getvalue('Refractory',ChAT);
Refractory_pChAT = getvalue('Refractory',pChAT);
Refractory_allChAT = getvalue('Refractory',allChAT);
BurstIndex_ChAT = getvalue('Burstindex',ChAT);
BurstIndex_pChAT = getvalue('Burstindex',pChAT);
BurstIndex_allChAT = getvalue('Burstindex',allChAT);
ThetaIndex_ChAT = getvalue('Thetaindex',ChAT);
ThetaIndex_pChAT = getvalue('Thetaindex',pChAT);
ThetaIndex_allChAT = getvalue('Thetaindex',allChAT);

% Tonic and phasic cells
phasic_ChAT = ChAT(Refractory_ChAT<=20);   % refractory below 20 ms
tonic_ChAT = ChAT(Refractory_ChAT>20);   % refractory above 20 ms
phasicB_ChAT = ChAT(Refractory_ChAT<=20&BurstIndex_ChAT>0.2);   % phasic, bursting
phasicNB_ChAT = ChAT(Refractory_ChAT<=20&BurstIndex_ChAT<=0.2);   % phasic, non-bursting
phasic_pChAT = pChAT(Refractory_pChAT<=20);   % refractory below 20 ms
tonic_pChAT = pChAT(Refractory_pChAT>20);   % refractory above 20 ms
phasicB_pChAT = pChAT(Refractory_pChAT<=20&BurstIndex_pChAT>0.2);   % phasic, bursting
phasicNB_pChAT = pChAT(Refractory_pChAT<=20&BurstIndex_pChAT<=0.2);   % phasic, non-bursting
phasic_allChAT = allChAT(Refractory_allChAT<=20);   % refractory below 20 ms
tonic_allChAT = allChAT(Refractory_allChAT>20);   % refractory above 20 ms
phasicB_allChAT = allChAT(Refractory_allChAT<=20&BurstIndex_allChAT>0.2);   % phasic, bursting
phasicNB_allChAT = allChAT(Refractory_allChAT<=20&BurstIndex_allChAT<=0.2);   % phasic, non-bursting

% Tonic and phasic cell indeces
phasic_ChATinx = Refractory_ChAT <= 20;   % refractory below 20 ms
tonic_ChATinx = Refractory_ChAT > 20;   % refractory above 20 ms
phasicB_ChATinx = Refractory_ChAT <= 20 & BurstIndex_ChAT > 0.2;   % phasic, bursting
phasicNB_ChATinx = Refractory_ChAT <= 20 & BurstIndex_ChAT <= 0.2;   % phasic, non-bursting
phasic_pChATinx = Refractory_pChAT <= 20;   % refractory below 20 ms
tonic_pChATinx = Refractory_pChAT > 20;   % refractory above 20 ms
phasicB_pChATinx = Refractory_pChAT <= 20 & BurstIndex_pChAT > 0.2;   % phasic, bursting
phasicNB_pChATinx = Refractory_pChAT <= 20 & BurstIndex_pChAT <= 0.2;   % phasic, non-bursting
phasic_allChATinx = Refractory_allChAT <= 20;   % refractory below 20 ms
tonic_allChATinx = Refractory_allChAT > 20;   % refractory above 20 ms
phasicB_allChATinx = Refractory_allChAT <= 20 & BurstIndex_allChAT > 0.2;   % phasic, bursting
phasicNB_allChATinx = Refractory_allChAT <= 20 & BurstIndex_allChAT <= 0.2;   % phasic, non-bursting


%%

figure
hold on
errorshade(1:32,mean(Spike_allChAT(phasic_allChATinx,:)),...
    nanse(Spike_allChAT(phasic_allChATinx,:)),...
    'LineColor',[0 0.8 0],'ShadeColor',[0 0.8 0])
errorshade(1:32,mean(Spike_NT),nanse(Spike_NT),'LineColor',[0 0 0],...
    'ShadeColor',[0 0 0])

%%

figure
hold on
errorshade(1:32,mean(Spike_ChAT(tonic_ChATinx,:)),...
    nanse(Spike_ChAT(tonic_ChATinx,:)),...
    'LineColor',green,'ShadeColor',green)
errorshade(1:32,mean(Spike_pChAT(tonic_pChATinx,:)),...
    nanse(Spike_pChAT(tonic_pChATinx,:)),...
    'LineColor',blue,'ShadeColor',blue)
errorshade(1:32,mean(Spike_NT),nanse(Spike_NT),'LineColor',[0 0 0],...
    'ShadeColor',[0 0 0])

%%

figure
hold on
errorshade(1:32,mean(Spike_ChAT(phasicNB_ChATinx,:)),...
    nanse(Spike_ChAT(phasicNB_ChATinx,:)),...
    'LineColor',green,'ShadeColor',green)
errorshade(1:32,mean(Spike_pChAT(phasicNB_pChATinx,:)),...
    nanse(Spike_pChAT(phasicNB_pChATinx,:)),...
    'LineColor',blue,'ShadeColor',blue)
errorshade(1:32,mean(Spike_NT),nanse(Spike_NT),'LineColor',[0 0 0],...
    'ShadeColor',[0 0 0])

%%

figure
hold on
errorshade(1:32,mean(Spike_ChAT(phasicB_ChATinx,:)),...
    nanse(Spike_ChAT(phasicB_ChATinx,:)),...
    'LineColor',green,'ShadeColor',green)
errorshade(1:32,mean(Spike_pChAT(phasicB_pChATinx,:)),...
    nanse(Spike_pChAT(phasicB_pChATinx,:)),...
    'LineColor',blue,'ShadeColor',blue)
errorshade(1:32,mean(Spike_NT),nanse(Spike_NT),'LineColor',[0 0 0],...
    'ShadeColor',[0 0 0])

%% other features

% Scatter plot
feature1 = 'PeakIntegralToIntegral';
feature2 = 'HalfPeakIntegralToIntegral';
figure
plot([SpikeShape_NT.(feature1)],[SpikeShape_NT.(feature2)],'o','MarkerSize',2,...
    'MarkerEdgeColor',grey,'MarkerFaceColor',grey)
hold on
plot([SpikeShape_ChAT.(feature1)],[SpikeShape_ChAT.(feature2)],'o',...
    'MarkerEdgeColor',green,'MarkerFaceColor',green)
plot([SpikeShape_pChAT.(feature1)],[SpikeShape_pChAT.(feature2)],'o',...
    'MarkerEdgeColor',blue,'MarkerFaceColor',blue)
xlabel(feature1)
ylabel(feature2)

% Stats
median([SpikeShape_NT.(feature1)])
median([SpikeShape_ChAT.(feature1)])
boxstat([SpikeShape_NT.(feature1)],[SpikeShape_ChAT.(feature1)],'NT','ChAT')
title(feature1)

median([SpikeShape_NT.(feature2)])
median([SpikeShape_ChAT.(feature2)])
boxstat([SpikeShape_NT.(feature2)],[SpikeShape_ChAT.(feature2)],'NT','ChAT')
title(feature2)

%% spike width vs depth

depth_ChAT = getvalue('DVpos',ChAT);
depth_pChAT = getvalue('DVpos',pChAT);
depth_allChAT = getvalue('DVpos',allChAT);
feature1 = 'PostValleyToPeak';

figure
hold on
plot(depth_ChAT,[SpikeShape_ChAT.(feature1)],'o',...
    'MarkerEdgeColor',green,'MarkerFaceColor',green)
plot(depth_pChAT,[SpikeShape_pChAT.(feature1)],'o',...
    'MarkerEdgeColor',blue,'MarkerFaceColor',blue)
xlabel('Depth')
ylabel(feature1)

% Regression
y = [SpikeShape_ChAT.(feature1)];
x = depth_ChAT';
[b,bint,r,rint,stats] = regress(y',[ones(length(x),1),x']);
R = sqrt(stats(1));         % correlation coefficient (R-value of the regression)
F = stats(2);           % F-test for H0: all coeff.-s are zero
p = stats(3);           % F-test significance
disp(p)

%% depth vs burstiness

depth_tonic_ChAT = getvalue('DVpos',tonic_ChAT);
depth_tonic_pChAT = getvalue('DVpos',tonic_pChAT);
depth_tonic_allChAT = getvalue('DVpos',tonic_allChAT);
depth_phasic_ChAT = getvalue('DVpos',phasic_ChAT);
depth_phasic_pChAT = getvalue('DVpos',phasic_pChAT);
depth_phasic_allChAT = getvalue('DVpos',phasic_allChAT);
depth_phasicB_ChAT = getvalue('DVpos',phasicB_ChAT);
depth_phasicB_pChAT = getvalue('DVpos',phasicB_pChAT);
depth_phasicB_allChAT = getvalue('DVpos',phasicB_allChAT);
depth_phasicNB_ChAT = getvalue('DVpos',phasicNB_ChAT);
depth_phasicNB_pChAT = getvalue('DVpos',phasicNB_pChAT);
depth_phasicNB_allChAT = getvalue('DVpos',phasicNB_allChAT);
feature1 = 'PostValleyToPeak';

figure
hold on
plot(depth_tonic_allChAT,[SpikeShape_allChAT(tonic_allChATinx).(feature1)],'o',...
    'MarkerEdgeColor',green,'MarkerFaceColor',green)
plot(depth_phasicB_allChAT,[SpikeShape_allChAT(phasicB_allChATinx).(feature1)],'o',...
    'MarkerEdgeColor',[0 0.8 0.8],'MarkerFaceColor',[0 0.8 0.8])
plot(depth_phasicNB_allChAT,[SpikeShape_allChAT(phasicNB_allChATinx).(feature1)],'o',...
    'MarkerEdgeColor',[0 0.2 0.8],'MarkerFaceColor',[0 0.2 0.8])
xlabel('Depth')
ylabel(feature1)

%% hit response

% Reward related firing rate change
Hitstat = getvalue('Hit_psth_stats',allChAT);   % hit trial PSTH statistics
HitPSTH = getvalue('Hit_psth',allChAT);   % hit trial PSTHs
[Hiteff_allChAT Hiteff2_allChAT baseline_allChAT] = deal(nan(1,length(allChAT)));
for k = 1:length(allChAT)
%     disp(allChAT{k})
%     pause(2)
    Hiteff_allChAT(k) = Hitstat{k}.maxvalue;   % maximal FR
    baseline_allChAT(k) = Hitstat{k}.baseline;  % baseline FR
    Hiteff2_allChAT(k) = log(Hitstat{k}.maxvalue/Hitstat{k}.baseline);  % relative FR change on log scale
end

% Reward related firing rate change
Hitstat = getvalue('Hit_psth_stats',pChAT);   % hit trial PSTH statistics
HitPSTH = getvalue('Hit_psth',pChAT);   % hit trial PSTHs
[Hiteff_pChAT Hiteff2_pChAT baseline_pChAT] = deal(nan(1,length(pChAT)));
for k = 1:length(pChAT)
    Hiteff_pChAT(k) = Hitstat{k}.maxvalue;   % maximal FR
    baseline_pChAT(k) = Hitstat{k}.baseline;  % baseline FR
    Hiteff2_pChAT(k) = log(Hitstat{k}.maxvalue/Hitstat{k}.baseline);  % relative FR change on log scale
end

% Reward related firing rate change
Hitstat = getvalue('Hit_psth_stats',ChAT);   % hit trial PSTH statistics
HitPSTH = getvalue('Hit_psth',ChAT);   % hit trial PSTHs
[Hiteff_ChAT Hiteff2_ChAT baseline_ChAT] = deal(nan(1,length(ChAT)));
for k = 1:length(ChAT)
    Hiteff_ChAT(k) = Hitstat{k}.maxvalue;   % maximal FR
    baseline_ChAT(k) = Hitstat{k}.baseline;  % baseline FR
    Hiteff2_ChAT(k) = log(Hitstat{k}.maxvalue/Hitstat{k}.baseline);  % relative FR change on log scale
end

%% depth vs Burstiness/spikewidth/postvalley

feature1 = 'PostValleyToPeak';
figure
hold on
plot(depth_ChAT,[SpikeShape_ChAT.(feature1)],'o',...
    'MarkerEdgeColor',green,'MarkerFaceColor',green)
plot(depth_pChAT,[SpikeShape_pChAT.(feature1)],'o',...
    'MarkerEdgeColor',blue,'MarkerFaceColor',blue)
xlabel('Depth')
ylabel(feature1)

feature1 = 'PeakToPostValleyTime';
figure
hold on
plot(depth_ChAT,[SpikeShape_ChAT.(feature1)],'o',...
    'MarkerEdgeColor',green,'MarkerFaceColor',green)
plot(depth_pChAT,[SpikeShape_pChAT.(feature1)],'o',...
    'MarkerEdgeColor',blue,'MarkerFaceColor',blue)
xlabel('Depth')
ylabel(feature1)

figure
hold on
plot(depth_ChAT,BurstIndex_ChAT,'o',...
    'MarkerEdgeColor',green,'MarkerFaceColor',green)
plot(depth_pChAT,BurstIndex_pChAT,'o',...
    'MarkerEdgeColor',blue,'MarkerFaceColor',blue)
xlabel('Depth')
ylabel('BurstIndex')

figure
hold on
plot(depth_ChAT,Refractory_ChAT,'o',...
    'MarkerEdgeColor',green,'MarkerFaceColor',green)
plot(depth_pChAT,Refractory_pChAT,'o',...
    'MarkerEdgeColor',blue,'MarkerFaceColor',blue)
xlabel('Depth')
ylabel('Refractory')

figure
hold on
plot(depth_ChAT,ThetaIndex_ChAT,'o',...
    'MarkerEdgeColor',green,'MarkerFaceColor',green)
plot(depth_pChAT,ThetaIndex_pChAT,'o',...
    'MarkerEdgeColor',blue,'MarkerFaceColor',blue)
xlabel('Depth')
ylabel('ThetaIndex')

feature1 = 'PostValleyToPeak';
figure
hold on
plot(BurstIndex_ChAT,[SpikeShape_ChAT.(feature1)],'o',...
    'MarkerEdgeColor',green,'MarkerFaceColor',green)
plot(BurstIndex_pChAT,[SpikeShape_pChAT.(feature1)],'o',...
    'MarkerEdgeColor',blue,'MarkerFaceColor',blue)
xlabel('Burst Index')
ylabel(feature1)

feature1 = 'PostValleyToPeak';
figure
hold on
plot(Refractory_ChAT,[SpikeShape_ChAT.(feature1)],'o',...
    'MarkerEdgeColor',green,'MarkerFaceColor',green)
plot(Refractory_pChAT,[SpikeShape_pChAT.(feature1)],'o',...
    'MarkerEdgeColor',blue,'MarkerFaceColor',blue)
xlabel('Refractory')
ylabel(feature1)

%% bursty vs non-bursty

feature1 = 'PostValleyToPeak';
median([SpikeShape_allChAT(tonic_allChATinx|phasicNB_allChATinx).(feature1)])
median([SpikeShape_allChAT(phasicB_allChATinx).(feature1)])
boxstat([SpikeShape_allChAT(tonic_allChATinx|phasicNB_allChATinx).(feature1)],[SpikeShape_allChAT(phasicB_allChATinx).(feature1)],'non-bursty','bursty')
title(feature1)
[H p] = ttest2([SpikeShape_allChAT(tonic_allChATinx|phasicNB_allChATinx).(feature1)],[SpikeShape_allChAT(phasicB_allChATinx).(feature1)])

median([depth_tonic_allChAT depth_phasicNB_allChAT])
median(depth_phasicB_allChAT)
boxstat([depth_tonic_allChAT' depth_phasicNB_allChAT'],depth_phasicB_allChAT,'non-bursty','bursty')
title('depth')
[H p] = ttest2([depth_tonic_allChAT' depth_phasicNB_allChAT'],depth_phasicB_allChAT)

figure
hold on
errorshade(1:32,mean(Spike_allChAT(tonic_allChATinx|phasicNB_allChATinx,:)),...
    nanse(Spike_allChAT(tonic_allChATinx|phasicNB_allChATinx,:)),...
    'LineColor',green,'ShadeColor',green)
errorshade(1:32,mean(Spike_allChAT(phasicB_allChATinx,:)),...
    nanse(Spike_allChAT(phasicB_allChATinx,:)),...
    'LineColor',[0 0.8 0.8],'ShadeColor',[0 0.8 0.8])
errorshade(1:32,mean(Spike_NT),nanse(Spike_NT),'LineColor',[0 0 0],...
    'ShadeColor',[0 0 0])

for k = 1:32
    boxstat(Spike_allChAT(tonic_allChATinx|phasicNB_allChATinx,k),Spike_allChAT(phasicB_allChATinx,k),'non-bursty','bursty')
end

%% bursty vs non-bursty #2

feature1 = 'PostValleyToPeak';
median([SpikeShape_allChAT(tonic_allChATinx).(feature1)])
median([SpikeShape_allChAT(phasicB_allChATinx).(feature1)])
boxstat([SpikeShape_allChAT(tonic_allChATinx).(feature1)],[SpikeShape_allChAT(phasicB_allChATinx).(feature1)],'non-bursty','bursty')
title(feature1)

median([depth_tonic_allChAT])
median(depth_phasicB_allChAT)
boxstat([depth_tonic_allChAT],depth_phasicB_allChAT,'non-bursty','bursty')
title('depth')

figure
hold on
errorshade(1:32,mean(Spike_allChAT(tonic_allChATinx,:)),...
    nanse(Spike_allChAT(tonic_allChATinx,:)),...
    'LineColor',green,'ShadeColor',green)
errorshade(1:32,mean(Spike_allChAT(phasicB_allChATinx,:)),...
    nanse(Spike_allChAT(phasicB_allChATinx,:)),...
    'LineColor',[0 0.8 0.8],'ShadeColor',[0 0.8 0.8])
errorshade(1:32,mean(Spike_NT),nanse(Spike_NT),'LineColor',[0 0 0],...
    'ShadeColor',[0 0 0])

%% depth vs stim-evoked FR (called baseline, but not really baseline)

% Punishment related firing rate change
FAstat = getvalue('FA_psth_stats',pChAT);   % hit trial PSTH statistics
FAPSTH = getvalue('FA_psth',pChAT);   % hit trial PSTHs
[FAeff_pChAT FAeff2_pChAT baseline_pChAT] = deal(nan(1,length(pChAT)));
for k = 1:length(pChAT)
    FAeff_pChAT(k) = FAstat{k}.maxvalue;   % maximal FR
    baseline_pChAT(k) = FAstat{k}.baseline;  % baseline FR
    FAeff2_pChAT(k) = log(FAstat{k}.maxvalue/FAstat{k}.baseline);  % relative FR change on log scale
end

% Punishment related firing rate change
FAstat = getvalue('FA_psth_stats',ChAT);   % hit trial PSTH statistics
FAPSTH = getvalue('FA_psth',ChAT);   % hit trial PSTHs
[FAeff_ChAT FAeff2_ChAT baseline_ChAT] = deal(nan(1,length(ChAT)));
for k = 1:length(ChAT)
    FAeff_ChAT(k) = FAstat{k}.maxvalue;   % maximal FR
    baseline_ChAT(k) = FAstat{k}.baseline;  % baseline FR
    FAeff2_ChAT(k) = log(FAstat{k}.maxvalue/FAstat{k}.baseline);  % relative FR change on log scale
end

figure
hold on
plot(depth_ChAT,baseline_ChAT,'o',...
    'MarkerEdgeColor',green,'MarkerFaceColor',green)
plot(depth_pChAT,baseline_pChAT,'o',...
    'MarkerEdgeColor',blue,'MarkerFaceColor',blue)
xlabel('Depth')
ylabel(feature1)

%% baseline FR vs SpikeWidth

Baseline_NT = getvalue('Baseline',NT);
Baseline_ChAT = getvalue('Baseline',ChAT);
Baseline_pChAT = getvalue('Baseline',pChAT);
Baseline_allChAT = getvalue('Baseline',allChAT);

%% baseline vs spike width

% Scatter plot
feature1 = 'PeakToPostValleyTime';
figure
plot([SpikeShape_NT.(feature1)]+rand(1,length(NT))*100-50,Baseline_NT,'o','MarkerSize',2,...
    'MarkerEdgeColor',grey,'MarkerFaceColor',grey)
hold on
plot([SpikeShape_ChAT.(feature1)],Baseline_ChAT,'o',...
    'MarkerEdgeColor',green,'MarkerFaceColor',green)
plot([SpikeShape_pChAT.(feature1)],Baseline_pChAT,'o',...
    'MarkerEdgeColor',blue,'MarkerFaceColor',blue)
xlabel(feature1)
ylabel('Firing Rate (Hz)')

% Stats
median([SpikeShape_NT.(feature1)])
median([SpikeShape_ChAT.(feature1)])
boxstat([SpikeShape_NT.(feature1)],[SpikeShape_ChAT.(feature1)],'NT','ChAT')
title(feature1)

median(Baseline_NT)
median(Baseline_ChAT)
boxstat(Baseline_NT,Baseline_ChAT,'NT','ChAT')
title('Firing Rate (Hz)')