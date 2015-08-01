function somresponseccg
%SOMRESPONSECCG   Response modulation strength summary.
%   SOMRESPONSECCG calculates average modulation strength for PSTHs around
%   behaviorally relevant events (cross-correlation between spike trains
%   and event times, CCG). Modulation strength is measured by SD of the
%   CCGs. SD values are normalized for each neuron to provide values
%   summing to 1 over the events. Average of significantly modulated CCGs
%   is taken for each event in the PV, SOM and non-tagged group.
%   Significance is assessed by crossing 0.01 significance levels
%   calculated based on shifted crosscorrelation (shifts ranged from 10 to
%   30s, see SOMRESPONSEPROFILES). Window size for CCGs was set to +-0.5
%   s and a bin width of 50 ms was used.
%
%   See also SOMRESPONSEPROFILES.

% Load
global DATADIR
global DATAPATH
load([DATADIR 'SOM_Sachin\PSTH\tagged_list_agreed'],'som_beh','pv_beh');
loadcb
load([DATAPATH 'SOM\ResponseCCG3\CCG_STRENGTH.mat'])

% Subsets
[is som_inx inx2] = intersect(CELLIDLIST,som_beh);   % indeces for SOM
[is pv_inx inx2] = intersect(CELLIDLIST,pv_beh);     % indeces for PV
nt_inx = setdiff((1:length(CELLIDLIST)),[som_inx pv_inx]);   % indeces for non-tagged

% Normalize and exclude non-significants
UaS = squeeze(CcgStrength(:,:,2));
LaS = squeeze(CcgStrength(:,:,3));
UaStrength = squeeze(CcgStrength(:,:,4));
for k = 1:size(CcgStrength,1)
    tuas = UaStrength(k,:);   % modulation strength
    tu = UaS(k,:);   % area of CCG above upper significance limit
    tl = LaS(k,:);   % area of CCG below lower significance limit
    tinx = tu<0.0001&tl<0.0001;   % indeces for non-significants (0.0001 is applied instead of 0 to deal with numerical errors)
    tuas(tinx) = NaN;   % non-significants set to NaN
    UaStrength(k,:) = tuas / nansum(tuas);   % normalize
end

% Plot groups
UaStrength1 = UaStrength(som_inx,:);    % SOM
figure;
plot(nansum(UaStrength1)/nansum(UaStrength1(:)),'Color','b')

UaStrength1 = UaStrength(pv_inx,:);     % PV
hold on
plot(nansum(UaStrength1)/nansum(UaStrength1(:)),'Color','r')

UaStrength1 = UaStrength(nt_inx,:);     % non-tagged
plot(nansum(UaStrength1)/nansum(UaStrength1(:)),'Color',[0.7 0.7 0.7])