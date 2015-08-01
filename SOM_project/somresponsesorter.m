function ppv_CellIDs = somresponsesorter
%SOMRESPONSESORTER   Find putative PV cells.
%   SOMRESPONSESORTER looks for putative PV cells in the 'non-port dataset'
%   based on response similarity, spike width and firing rate. Response
%   similarity to home-zone-out event is evaluated based on dot product
%   (equivalent to correlation) between smoothed peri-stimulus time
%   histograms.
%
%   PPV_CELLIDS = SOMRESPONSESORTER returns the cell IDs for putative
%   PV cells.
%
%   See also SOMRESPONSESORTER_PORT.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   6-Oct-2012

%   Edit log: BH 10/12/12

% Load HomeZoneOut data
global DATADIR
load([DATADIR 'SOM_Sachin\PSTH\HomeZout.mat'])
load([DATADIR 'SOM_Sachin\PSTH\SpikeWidth.mat'])
pv_Psth = Pv_psth;   % response profile (PSTH)
nt_Psth = Non_tagged_psth;
pv_SpikeWidth = Pv_SpikeWidth;  % spike width
nt_SpikeWidth = Non_tagged_SpikeWidth;
pv_FiringRate = Pv_Rate;  % firing rate
nt_FiringRate = Non_tagged_Rate;
pv_CellIDs = Pv_cells;   % cell IDs
nt_CellIDs = Non_tagged_cells;

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

mspvd = mean(spvd1);  % mean similarity for PVs
ntspvdd = abs(sntd1-mspvd);   % distance from PV mean
pvspvdd = abs(spvd1-mspvd);
figure   % similarity distance from mean PV vs firing rate
plot(ntspvdd,nt_FiringRate,'o','Color',[0.7 0.7 0.7])
hold on
plot(pvspvdd,pv_FiringRate,'ro','MarkerFaceColor','red','MarkerEdgeColor','red')
xlabel('Correlation with mean PV PSTH')
ylabel('Firing rate')

% Distribution of similarity
figure
hist([sntd1(nt_SpikeWidth<270) spvd1(PvSpikeWidth<270)],25)   % restrict to narrow spiking cells
P1 = findobj(gca,'type','patch');
set(P1,'FaceColor','w')
[nm xout] = hist([sntd1(nt_SpikeWidth<270) spvd1(pv_SpikeWidth<270)],25);  % get binning information
hold on
hist(spvd1(pv_SpikeWidth<270),xout)   % overlay PVs similarity
P = findobj(gca,'type','patch');
P2 = setdiff(P,P1);
set(P2,'FaceColor','r')
nm = hist(sntd1(nt_SpikeWidth>270),xout);   % overlay wide spiking cells
plot(xout,nm/max(nm)*10,'Color',[0.7 0.7 0.7])

figure
hist([sntd1(nt_FiringRate>15&nt_SpikeWidth<270) spvd1(pv_FiringRate>15&pv_SpikeWidth<270)],18)   % restrict to narrow spiking fast firing cells
P1 = findobj(gca,'type','patch');
set(P1,'FaceColor','w')
[nm xout] = hist([sntd1(nt_FiringRate>15&nt_SpikeWidth<270) spvd1(pv_FiringRate>15&pv_SpikeWidth<270)],18);  % get binning information
hold on
hist(spvd1(pv_FiringRate>15&pv_SpikeWidth<270),xout)   % overlay PVs similarity
P = findobj(gca,'type','patch');
P2 = setdiff(P,P1);
set(P2,'FaceColor','r')
nm = hist(sntd1(nt_FiringRate<15|nt_SpikeWidth>270),xout);   % overlay wide spiking/slow firing cells
plot(xout,nm/max(nm)*10,'Color',[0.7 0.7 0.7])

% Plot non-tagged cells sorted by similarity
[srt ia] = sort(sntd1,'descend');  % sort
figure
imagesc(snt1(ia,:))

% Find putative PV cells
pPVinx = sntd1'>60&nt_SpikeWidth<270&nt_FiringRate>15;  % use similarity and spike width cutoff
ppv_CellIDs = nt_CellIDs(pPVinx);

keyboard

figure
hist([sntd1(nt_FiringRate>12) spvd1(pv_FiringRate>12)],25)   % restrict to narrow spiking cells
P1 = findobj(gca,'type','patch');
set(P1,'FaceColor','w')
[nm xout] = hist([sntd1(nt_FiringRate>12) spvd1(pv_FiringRate>12)],25);  % get binning information
hold on
hist(spvd1(pv_FiringRate>12),xout)   % overlay PVs similarity
P = findobj(gca,'type','patch');
P2 = setdiff(P,P1);
set(P2,'FaceColor','r')

figure
hist([ntspvdd(nt_FiringRate>12) pvspvdd(pv_FiringRate>12)],25)   % restrict to narrow spiking cells
P1 = findobj(gca,'type','patch');
set(P1,'FaceColor','w')
[nm xout] = hist([ntspvdd(nt_FiringRate>12) pvspvdd(pv_FiringRate>12)],25);  % get binning information
hold on
hist(pvspvdd(pv_FiringRate>12),xout)   % overlay PVs similarity

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