function nbtonicphasic
%NBTONICPHASIC   Group cells based on auto-correlation.
%   Cholinergic and putative cholinergic cells are grouped into tonic,
%   phasic - bursting and phasic - non-bursting groups based on refractory
%   and burst index calculated from their auto-correlograms (see NBACG for
%   details). ACGs are smoothed and normaized before averaging.
%
%   See also NBACG.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   14-June-2013

% Load auto-correlations
load('C:\Balazs\_analysis\NB\ACG\medium\ACG_matrices.mat')

% Color code
allChAT = cellids;
NumChAT = length(allChAT);   % number of cells
mouse = getvalue('RatID',allChAT);  % animal ID
nms = length(unique(mouse));   % number of mice (n=7)
mouse2 = nan(size(mouse));
for k = 1:nms
    mouse2(mouse==min(mouse)) = k;   % label for indiv. mice
    mouse(mouse==min(mouse)) = Inf;  % temp. variable
end
clr = colormap(jet(nms));   % unique color for each mouse

ChAT = getvalue('ChAT+',allChAT);
pChAT = getvalue('pChAT+',allChAT);
clr2 = [ChAT>0 pChAT>0 zeros(NumChAT,1)];   % coloring according to identified and putative ChATs

% Plot ACG properties
H1 = figure;
hold on
xlabel('Refractory (ms)')
ylabel('Burst index')
H2 = figure;
hold on
xlabel('Refractory (ms)')
ylabel('Burst index')
for k = 1:NumChAT
    figure(H1)
    plot(Refractory(k),BurstIndex(k),'o','MarkerSize',12,'MarkerFaceColor',clr(mouse2(k),:),'MarkerEdgeColor','k')  % scatter plot
    figure(H2)
    plot(Refractory(k),BurstIndex(k),'o','MarkerSize',12,'MarkerFaceColor',clr2(k,:),'MarkerEdgeColor','k')  % scatter plot
end

% Normalize ACG
CCRnorm = nan(size(CCR));
for k = 1:NumChAT
    sccr = smooth(CCR(k,:),'linear',5);   % smooth
%     CCRnorm(k,:) = zscore(sccr);   % z-score
    b = mean(sccr(lags>-500&lags<500));
    CCRnorm(k,:) = sccr / b;
end

% Tonic and phasic cells
phasic = Refractory <= 20;   % refractory below 20 ms
tonic = Refractory > 20;   % refractory above 20 ms
phasicB = Refractory <= 20 & BurstIndex > 0.2;   % phasic, bursting
phasicNB = Refractory <= 20 & BurstIndex <= 0.2;   % phasic, non-bursting

% Plot tonice cells
tonicCCR = CCRnorm(tonic,:);  % ACG for tonic cells
[srt inx] = sort(Refractory(tonic),'descend');   % sort according to refractory
figure
imagesc(lags(lags>-250&lags<250),1:size(tonicCCR,1),tonicCCR(inx,lags>-250&lags<250))
set(gca,'CLim',[0 2.5])
title('ACG of TONIC cholinergic cells, sorted by Refractory')

% Plot phasic cells
phasicBCCR = CCRnorm(phasicB,:);  % ACG for phasic bursting cells
[srt inx] = sort(BurstIndex(phasicB),'descend');   % sort according to burst index
figure
imagesc(lags(lags>-250&lags<250),1:size(phasicBCCR,1),phasicBCCR(inx,lags>-250&lags<250))
set(gca,'CLim',[0 2.5])
title('ACG of PHASIC BURSTING cholinergic cells, sorted by BurstIndex')

phasicNBCCR = CCRnorm(phasicNB,:);  % ACG for phasic non-bursting cells
[srt inx] = sort(BurstIndex(phasicNB),'descend');   % sort according to burst index
figure
imagesc(lags(lags>-250&lags<250),1:size(phasicNBCCR,1),phasicNBCCR(inx,lags>-250&lags<250))
set(gca,'CLim',[0 2.5])
title('ACG of PHASIC NON-BURSTING cholinergic cells, sorted by BurstIndex')

% Plot averages
mn_tonic = mean(tonicCCR(:,lags>-250&lags<250));   % average ACG, tonic cells
se_tonic = std(tonicCCR(:,lags>-250&lags<250)) / sqrt(size(tonicCCR,1));   % SE, tonic cells
mn_phasicB = mean(phasicBCCR(:,lags>-250&lags<250));   % average ACG, phasic, bursting cells
se_phasicB = std(phasicBCCR(:,lags>-250&lags<250)) / sqrt(size(phasicB,1));   % SE, phasic, bursting cells
mn_phasicNB = mean(phasicNBCCR(:,lags>-250&lags<250));   % average ACG, phasic, non-bursting cells
se_phasicNB = std(phasicNBCCR(:,lags>-250&lags<250)) / sqrt(size(phasicNB,1));   % SE, non-bursting cells

figure
hold on
errorshade(lags(lags>-250&lags<250),mn_tonic,se_tonic,'LineColor',[0 0.8 0],'ShadeColor',[0 0.8 0],'LineWidth',2)
errorshade(lags(lags>-250&lags<250),mn_phasicB,se_phasicB,'LineColor',[0 0.8 0.8],'ShadeColor',[0 0.8 0.8],'LineWidth',2)
errorshade(lags(lags>-250&lags<250),mn_phasicNB,se_phasicNB,'LineColor',[0 0.2 0.8],'ShadeColor',[0 0.2 0.8],'LineWidth',2)
xlim([-250 250])
xlabel('Time (ms)')
ylabel('Normalized count')

keyboard