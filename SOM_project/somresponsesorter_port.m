function ppv_CellIDs = somresponsesorter_port
%SOMRESPONSESORTER_PORT   Find putative PV cells.
%   SOMRESPONSESORTER_PORT looks for putative PV cells in the 'port
%   dataset' based on response similarity, spike width and firing rate.
%   Response similarity to port-out event is evaluated based on dot product
%   (equivalent to correlation) between smoothed peri-stimulus time
%   histograms.
%
%   PPV_CELLIDS = SOMRESPONSESORTER_PORT returns the cell IDs for putative
%   PV cells.
%
%   See also SOMRESPONSESORTER.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   6-Oct-2012

%   Edit log: BH 10/11/12

% Load PortOut data
global DATADIR
load([DATADIR 'SOM_Sachin\PSTH_port\PortOut.mat'])
pv_Psth = pv;   % response profile (PSTH)
nt_Psth = non_tagged;
pv_SpikeWidth = Pv_Spikewidth;  % spike width
nt_SpikeWidth = NonTaggedAllSpikeWidth;
pv_FiringRate = Pv_Rate;  % firing rate
nt_FiringRate = Rate;
pv_CellIDs = pv_cellids;   % cell IDs
nt_CellIDs = non_tagged_cellids;

% Remove NaNs
[x1 y1] = find(isnan(nt_Psth));  % NaN rows in Psth
x2 = find(isnan(nt_SpikeWidth));  % NaNs in spike width vector
x = union(unique(x1),unique(x2));
nt_Psth(x,:) = [];
nt_SpikeWidth(x) = [];
nt_FiringRate(x) = [];
nt_CellIDs(x) = [];

% Restrict in time
pv_Psth = pv_Psth(:,7:207);
nt_Psth = nt_Psth(:,7:207);

% Smooth
spv1 = smoothallcells(pv_Psth);   % PV, PortOut
snt1 = smoothallcells(nt_Psth);   % non-tagged, PortOut

% Correlation with specific PSTHs
template = mean(spv1);   % correlation with mean PV, PortOut
spvd1 = tcorr(spv1,template);    % PV, PortOut
sntd1 = tcorr(snt1,template);    % non-tagged, PortOut

% Spike width distribution
figure
hist(nt_SpikeWidth,65)

% Scatter plots for similarity, spike width and firing rate
figure   % spike width vs similarity
plot(sntd1,nt_SpikeWidth,'o','Color',[0.7 0.7 0.7])
hold on
plot(spvd1,pv_SpikeWidth,'ro','MarkerFaceColor','red','MarkerEdgeColor','red')
xlabel('Correlation with mean PV PSTH')
ylabel('Spike width')

mpvsw = mean(pv_SpikeWidth(pv_SpikeWidth<270));  % mean spike width of narrow spiking PVs
ntswd = abs(nt_SpikeWidth-mpvsw);   % distance from PV mean
pvswd = abs(pv_SpikeWidth-mpvsw);
figure   % spike width distance from mean PV vs similarity
plot(sntd1,ntswd,'o','Color',[0.7 0.7 0.7])
hold on
plot(spvd1,pvswd,'ro','MarkerFaceColor','red','MarkerEdgeColor','red')
xlabel('Correlation with mean PV PSTH')
ylabel('Spike width')

figure   % similarity vs. firing rate
plot(sntd1,nt_FiringRate,'o','Color',[0.7 0.7 0.7])
hold on
plot(spvd1,pv_FiringRate,'ro','MarkerFaceColor','red','MarkerEdgeColor','red')
xlabel('Correlation with mean PV PSTH')
ylabel('Firing rate')

% Distribution of similarity
figure
hist([sntd1(nt_SpikeWidth<270) spvd1(pv_SpikeWidth<270)],25)   % restrict to narrow spiking cells
P1 = findobj(gca,'type','patch');
set(P1,'FaceColor','w')
[nm xout] = hist([sntd1(nt_SpikeWidth<270) spvd1(pv_SpikeWidth<270)],25);  % get binning information
hold on
hist(spvd1(pv_SpikeWidth<270),xout)   % overlay PVs similarity

% Plot non-tagged cells sorted by similarity
[srt ia] = sort(sntd1,'descend');  % sort
figure
imagesc(snt1(ia,:))

% Find putative PV cells
similarity_limit = 90;
spike_width_limit = 270;
firing_rate_limit = 15;
pPVinx = sntd1'>similarity_limit&nt_SpikeWidth<spike_width_limit&nt_FiringRate'>firing_rate_limit;  % use similarity and spike width cutoff
ppv_CellIDs = nt_CellIDs(pPVinx);

keyboard

% Illustration for selection of putative PVs
figure   % similarity vs. firing rate
plot(sntd1,nt_FiringRate,'o','Color',[0.7 0.7 0.7])
hold on
plot(spvd1,pv_FiringRate,'o','MarkerFaceColor',[0.8 0 0],'MarkerEdgeColor',[0.8 0 0])
plot(sntd1(pPVinx),nt_FiringRate(pPVinx),'o','MarkerFaceColor',[1 0.6 0],'MarkerEdgeColor',[1 0.6 0])
xlabel('Correlation with mean PV PSTH')
ylabel('Firing rate')
line([similarity_limit similarity_limit],ylim,'LineStyle',':','Color','k')
line(xlim,[firing_rate_limit firing_rate_limit],'LineStyle',':','Color','k')
set(gca,'TickDir','out')
box off

figure   % spike width vs. firing rate
plot(nt_SpikeWidth,nt_FiringRate,'o','Color',[0.7 0.7 0.7])
hold on
plot(pv_SpikeWidth,pv_FiringRate,'o','MarkerFaceColor',[0.8 0 0],'MarkerEdgeColor',[0.8 0 0])
plot(nt_SpikeWidth(pPVinx),nt_FiringRate(pPVinx),'o','MarkerFaceColor',[1 0.6 0],'MarkerEdgeColor',[1 0.6 0])
xlabel('Spike width')
ylabel('Firing rate')
line([spike_width_limit spike_width_limit],ylim,'LineStyle',':','Color','k')
line(xlim,[firing_rate_limit firing_rate_limit],'LineStyle',':','Color','k')
set(gca,'TickDir','out')
box off

figure   % spike width vs. similarity
plot(nt_SpikeWidth,sntd1,'o','Color',[0.7 0.7 0.7])
hold on
plot(pv_SpikeWidth,spvd1,'o','MarkerFaceColor',[0.8 0 0],'MarkerEdgeColor',[0.8 0 0])
plot(nt_SpikeWidth(pPVinx),sntd1(pPVinx),'o','MarkerFaceColor',[1 0.6 0],'MarkerEdgeColor',[1 0.6 0])
xlabel('Spike width')
ylabel('Correlation with mean PV PSTH')
line([spike_width_limit spike_width_limit],ylim,'LineStyle',':','Color','k')
line(xlim,[similarity_limit similarity_limit],'LineStyle',':','Color','k')
set(gca,'TickDir','out')
box off

% -------------------------------------------------------------------------
function S = smoothallcells(X)

% Smoothing with moving average
sx = size(X);
S = zeros(sx);
for k = 1:sx(1)
    S(k,:) = smooth(X(k,:),'linear',11);
end

% -------------------------------------------------------------------------
function C = tcorr(X,template)

% Correlation with specific PSTHs
C = template * X';