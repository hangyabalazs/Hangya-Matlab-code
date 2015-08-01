function condPerf = cond_accuracy_fr2(lim1,lim2)
%COND_ACCURACY_FR2   Accuracy conditioned on firing rate.
%   COND_ACCURACY_FR2 calculates the percentage of hits relative to the sum
%   of hits and false alarms for specific firing rate intervals of a given
%   cell. Firing rates are binned; hits and false alarms are partitioned
%   according to the firing rate bins. It does not handle multiple
%   sessions; different stimulus intensities are pooled together. To a more
%   sophisticated, however slower solution, see COND_ACCURACY_FR.
%
%   CONDPERF = COND_ACCURACY_FR2(LIM1,LIM2) uses LIM1 and LIM2 to determine
%   the time window before stimuli used for firing rate calculation (in
%   seconds). Default is -2:0. Conditional performance vector (CONDPERF) is
%   returned.
%
%   See also COND_ACCURACY_FR.

%   Edit log: BH 7/7/11

% Input variables
if nargin == 0
    lim1 = -2;
    lim2 = 0;
end
animalID = 'n013';
sessionIDs = {'110628a'};
NUMsessions = length(sessionIDs);
sessionID = sessionIDs{1};
fullpth = [getpref('cellbase','datapath') '\' animalID '\' sessionID '\'];
TE = load([fullpth 'TrialEvents.mat']);
load([fullpth 'EVENTSPIKES5_4_long2.mat'])
spikes_stimon = event_stimes{1};

% Prestimulus frequency
NUMtrials = length(spikes_stimon);
prestimfreq = nan(1,NUMtrials);
for k = 1:NUMtrials
    lspikes = spikes_stimon{k};
    lspikes2 = lspikes(lspikes>lim1&lspikes<lim2);   % time window: one sec before stimulus onset
    prestimfreq(k) = length(lspikes2) / (lim2 - lim1);
end

% Response types
ishit = logical(nan2zero(TE.Hit));
ismiss = logical(nan2zero(TE.Miss));
isfa = logical(nan2zero(TE.FalseAlarm));
iscr = logical(nan2zero(TE.CorrectRejection));

% Conditional performance
edges = (10:10:40);
cnts = (edges(1:end-1) + edges(2:end)) / 2;
NUMbins = length(cnts);
condPerf = zeros(1,NUMbins);
for iB = 1:NUMbins
    binhits = prestimfreq(ishit) >= edges(iB) & prestimfreq(ishit) < edges(iB+1);
    nobh = sum(binhits);
    binfas = prestimfreq(isfa) >= edges(iB) & prestimfreq(isfa) < edges(iB+1);
    nobf = sum(binfas);
    
    condPerf(iB) = nobh / (nobh + nobf);    % hits / (hits + FAs)
end
% figure    % plot
plot(cnts,condPerf,'Color',[0 (60+lim1)/60 0])