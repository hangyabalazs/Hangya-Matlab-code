function [H C] = cond_accuracy_fr3_test(cellid,varargin)
%COND_ACCURACY_FR3   Accuracy conditioned on firing rate.
%   COND_ACCURACY_FR3(CELLID) calculates the percentage of correct
%   responses (accuracy) and the percentage of response trials (hits and/or 
%   false alarms) alarms for specific firing rate intervals of a given cell
%   (CELLID). Firing rates are binned; hits, false alarms, misses and 
%   correct rejections are partitioned according to the firing rate bins.
%   Different accuracy measures (go performance, no-go performance, overall
%   performance and hit/hit+FA rate, termed 'action performance') and
%   response rate measures (go respose rate, i.e go performance; no-go
%   response rate, overall respose rate, respose rate in correct trials)
%   are plotted agains reaction time.
%
%   Trials are partitioned to low (10-50 percentile) and high (50-90 
%   percentile) reaction time; performance is plotted conditioned on low
%   and high reaction time.
%
%   Trials are partitioned to low and high firing rate (median split).
%   Performance and response rate measures are plotted for low and high
%   firing rate.
%
%   [H CONDPERF] = COND_ACCURACY_FR3(CELLID) returns handles of the figures
%   (H) and all performance measures (CONDPERF).
%
%   Optional input parameter-value pairs (with default values):
%       'window', [0 0.5] - firing rate window relative to the reference
%           event
%       'event', 'LeftPortIn' - reference event for firing rate window
%       'limit2iti', false - limit the firing rate window to the foreperiod
%           if the window size is larger; for 'StimulusOn' reference event
%       'display', false - controls plotting
%
%   See also COND_ACCURACY_FR and COND_ACCURACY_FR2.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   29-Oct-2013

%   Edit log: BH 10/29/13

% Default arguments
prs = inputParser;
addRequired(prs,'cellid',@iscellid)
addParamValue(prs,'window',[0 0.5],@(s)isnumeric(s)&isequal(length(s),2))  % time window relative to the event, in seconds
addParamValue(prs,'event','LeftPortIn',@ischar)   % default reference event: 'LeftPortIn'
addParamValue(prs,'limit2iti',false,@(s)islogical(s)|ismember(s,[0 1]))   % restrict firing rate window to foreperiod
addParamValue(prs,'display',false,@(s)islogical(s)|ismember(s,[0 1]))   % control displaying rasters and PSTHs
parse(prs,cellid,varargin{:})
g = prs.Results;

% Load trial events
ST = loadcb(cellid,'EVENTSPIKES');   % load prealigned spikes for stimulation events
TE = loadcb(cellid,'TrialEvents');   % load trial events
event_pos = findcellstr(ST.events(:,1),g.event);  % select event
spikes_stimon = ST.event_stimes{event_pos};  % spike times

% Prestimulus frequency
NUMtrials = length(spikes_stimon);   % number of trials
prestimfreq = nan(1,NUMtrials);
lim1 = g.window(1);   % window boundaries
lim2 = g.window(2);
for k = 1:NUMtrials   % loop through trials
    if g.limit2iti
        itiwin = TE.ITIDistribution(k);   % ITI
        llim1 = max(-itiwin,lim1);
    else
        llim1 = lim1;
    end
    llim2 = lim2;
    lspikes = spikes_stimon{k};   % spike times in the current trial
    lspikes2 = lspikes(lspikes>llim1&lspikes<llim2);   % apply time window
    prestimfreq(k) = length(lspikes2) / (llim2 - llim1);   % firing rate
end

% Response types
ishit = logical(nan2zero([TE.Hit(1:end-1) 0]));   % Hit trials
ismiss = logical(nan2zero([TE.Miss(1:end-1) 0]));   % Miss trials
isfa = logical(nan2zero([TE.FalseAlarm(1:end-1) 0]));   % False alarm trials
iscr = logical(nan2zero([TE.CorrectRejection(1:end-1) 0]));   % Correct rejection trials

% Split psychometric function based on prestimulus FR
prcthres = [0 50];   % percentile threshold
frthres = prctile(prestimfreq,prcthres);   % FR threshold
lowFR = find(prestimfreq>=frthres(1)&prestimfreq<frthres(2));   % FRs of the percentile interval
numLFR = length(lowFR);   % number of low FR trials
prcthres = [50 100];   % percentile threshold
frthres = prctile(prestimfreq,prcthres);   % FR threshold
highFR = find(prestimfreq>frthres(1)&prestimfreq<=frthres(2));   % FRs of the percentile interval
numHFR = length(highFR);   % number of high FR trials
C = mainCondPerf(TE,lowFR,highFR,ishit,ismiss,isfa,iscr);

% Permutations
pno = 1000;   % number of permutations
allFR = [lowFR highFR];
for iP = 1:pno
    allFRt = allFR(randperm(numLFR+numHFR));
    lowFRt = allFRt(1:numLFR);   % randomized 'low FR' trials
    highFRt = allFRt(numLFR+1:end);   % randomized 'high FR' trials
    Cperm(iP) = mainCondPerf(TE,lowFRt,highFRt,ishit,ismiss,isfa,iscr);
end

% Test
origD = C.condDiscriminationS - C.condDiscriminationL;
permD1 = {Cperm.condDiscriminationS};
permD1 = cell2mat(permD1');
permD2 = {Cperm.condDiscriminationL};
permD2 = cell2mat(permD2');
permD = permD1 - permD2;
p_values = sum(permD-repmat(origD,pno,1)>0) / pno;

% Plot
if g.display
%     H.Hfra = figure('Color',[0.3 0.3 0.3],'Position',[456 314 840 664]);   % accuracy measures
%     subplot(221)    % plot (all trials)
%     clr = [0.5 0 0; 1 0 0; 1 0.5 0; 1 1 0];
%     plot(SI,C.condPerfS,'Color',clr(1,:),'LineWidth',3)
%     hold on
%     plot(SI,C.condPerfL,'Color',clr(4,:),'LineWidth',3)
%     set(gca,'Color',[0.8 0.8 0.8],'LineWidth',2)
%     title('all trials')
%     xlabel('SPL (dB)')
%     ylabel('% correct')
%     
%     subplot(222)    % plot (go trials)
%     clr = [0.5 0 0; 1 0 0; 1 0.5 0; 1 1 0];
%     plot(SI,C.condGoPerfS,'Color',clr(1,:),'LineWidth',3)
%     hold on
%     plot(SI,C.condGoPerfL,'Color',clr(4,:),'LineWidth',3)
%     set(gca,'Color',[0.8 0.8 0.8],'XColor',[0 0.8 0],'YColor',[0 0.8 0],'LineWidth',2)
%     title('go trials','Color',[0 0.8 0])
%     xlabel('SPL (dB)')
%     ylabel('% correct')
%     
%     subplot(223)    % plot (no-go trials)
%     clr = [0.5 0 0; 1 0 0; 1 0.5 0; 1 1 0];
%     plot(SI,C.condNoGoPerfS,'Color',clr(1,:),'LineWidth',3)
%     hold on
%     plot(SI,C.condNoGoPerfL,'Color',clr(4,:),'LineWidth',3)
%     set(gca,'Color',[0.8 0.8 0.8],'XColor',[0.8 0 0],'YColor',[0.8 0 0],'LineWidth',2)
%     title('no-go trials','Color',[0.8 0 0])
%     xlabel('SPL (dB)')
%     ylabel('% correct')
%     
%     subplot(224)    % plot (response trials)
%     clr = [0.5 0 0; 1 0 0; 1 0.5 0; 1 1 0];
%     plot(SI,C.condActionPerfS,'Color',clr(1,:),'LineWidth',3)
%     hold on
%     plot(SI,C.condActionPerfL,'Color',clr(4,:),'LineWidth',3)
%     set(gca,'Color',[0.8 0.8 0.8],'XColor',[0.8 0.8 0],'YColor',[0.8 0.8 0],'LineWidth',2)
%     title('response trials','Color',[0.8 0.8 0])
%     xlabel('SPL (dB)')
%     ylabel('% correct')
%     
%     H.Hfrr = figure('Color',[0.3 0.3 0.3],'Position',[456 314 840 664]);   % action (response) measures
%     subplot(221)    % plot (all trials)
%     clr = [0.5 0 0; 1 0 0; 1 0.5 0; 1 1 0];
%     plot(SI,C.condRespS,'Color',clr(1,:),'LineWidth',3)
%     hold on
%     plot(SI,C.condRespL,'Color',clr(4,:),'LineWidth',3)
%     set(gca,'Color',[0.8 0.8 0.8],'LineWidth',2)
%     title('all trials')
%     xlabel('SPL (dB)')
%     ylabel('% response')
%     
%     subplot(222)    % plot (go trials)
%     clr = [0.5 0 0; 1 0 0; 1 0.5 0; 1 1 0];
%     plot(SI,C.condGoRespS,'Color',clr(1,:),'LineWidth',3)
%     hold on
%     plot(SI,C.condGoRespL,'Color',clr(4,:),'LineWidth',3)
%     set(gca,'Color',[0.8 0.8 0.8],'XColor',[0 0.8 0],'YColor',[0 0.8 0],'LineWidth',2)
%     title('go trials','Color',[0 0.8 0])
%     xlabel('SPL (dB)')
%     ylabel('% response')
%     
%     subplot(223)    % plot (no-go trials)
%     clr = [0.5 0 0; 1 0 0; 1 0.5 0; 1 1 0];
%     plot(SI,C.condNoGoRespS,'Color',clr(1,:),'LineWidth',3)
%     hold on
%     plot(SI,C.condNoGoRespL,'Color',clr(4,:),'LineWidth',3)
%     set(gca,'Color',[0.8 0.8 0.8],'XColor',[0.8 0 0],'YColor',[0.8 0 0],'LineWidth',2)
%     title('no-go trials','Color',[0.8 0 0])
%     xlabel('SPL (dB)')
%     ylabel('% response')
%     
%     subplot(224)    % plot (correct trials)
%     clr = [0.5 0 0; 1 0 0; 1 0.5 0; 1 1 0];
%     plot(SI,C.condCorrectRespS,'Color',clr(1,:),'LineWidth',3)
%     hold on
%     plot(SI,C.condCorrectRespL,'Color',clr(4,:),'LineWidth',3)
%     set(gca,'Color',[0.8 0.8 0.8],'XColor',[0.8 0.8 0],'YColor',[0.8 0.8 0],'LineWidth',2)
%     title('correct trials','Color',[0.8 0.8 0])
%     xlabel('SPL (dB)')
%     ylabel('% response')
%     
    H.Hdpr = figure('Color',[0.3 0.3 0.3],'Position',[456 314 840 664]);   % action (response) measures
    subplot(211)    % plot (hit and false alarm rates)
    plot(SI,C.condGoRespS,'Color',clr(1,:),'LineWidth',3)
    hold on
    plot(SI,C.condGoRespL,'Color',clr(4,:),'LineWidth',3)
    plot(SI,C.condNoGoRespS,'Color',clr(1,:),'LineWidth',3,'LineStyle',':')
    plot(SI,C.condNoGoRespL,'Color',clr(4,:),'LineWidth',3,'LineStyle',':')
    set(gca,'LineWidth',2)
    legend({'Hit rate, high FR' 'Hit rate, low FR' 'FA rate, high FR' 'FA rate, low FR'},...
        'Location','best')
    title('Performance')
    xlabel('SPL (dB)')
    ylabel('% hit (solid) or FA (dashed) response')
    
    subplot(223)    % plot (hit-false alarm)
    plot(SI,C.condDiscriminationS-C.condDiscriminationL,'Color','k','LineWidth',3)
    set(gca,'LineWidth',2)
    title('hit rate - false alarm rate')
    xlabel('SPL (dB)')
    ylabel('% hit - FA response (high FR - low FR)')
    
    subplot(224)    % plot (d')
    plot(SI,C.condDPrimeS-C.condDPrimeL,'Color','k','LineWidth',3)
    set(gca,'LineWidth',2)
    title('d''')
    xlabel('SPL (dB)')
    ylabel('d'' (high FR - low FR)')
end

% -------------------------------------------------------------------------
function C = mainCondPerf(TE,lowFR,highFR,ishit,ismiss,isfa,iscr)

% Find unique sound intensities
SI = unique(TE.StimulusDuration);
SI = SI(~isnan(SI));
numSI = length(SI);   % number of different stimulus intensities
[C.condPerfS C.condGoPerfS C.condNoGoPerfS C.condActionPerfS ...
    C.condRespS C.condGoRespS C.condNoGoRespS C.condCorrectRespS ...
    C.condDiscriminationS C.condDPrimeS nobS nobhS nobfS nobmS nobcS ...
    C.condPerfL C.condGoPerfL C.condNoGoPerfL C.condActionPerfL ...
    C.condRespL C.condGoRespL C.condNoGoRespL C.condCorrectRespL ...
    C.condDiscriminationL C.condDPrimeL nobL nobhL nobfL nobmL nobcL] = ...
    deal(zeros(1,numSI));

for iB = 1:numSI   % loop through stim. intensities
    bins = TE.StimulusDuration == SI(iB);
    nobS(iB) = length(intersect(find(bins),highFR));   % number of trials in the bin
    binhits = intersect(find(ishit&bins),highFR);
    nobhS(iB) = length(binhits);   % number of Hit trials in the bin
    binfas = intersect(find(isfa&bins),highFR);
    nobfS(iB) = length(binfas);   % number of False alarm trials in the bin
    binmisss = intersect(find(ismiss&bins),highFR);
    nobmS(iB) = length(binmisss);   % number of Miss trials in the bin
    bincrs = intersect(find(iscr&bins),highFR);
    nobcS(iB) = length(bincrs);   % number of Correct rejection trials in the bin
    C.condPerfS(iB) = (nobhS(iB) + nobcS(iB)) / nobS(iB);    % percent correct
    C.condGoPerfS(iB) = nobhS(iB) / (nobhS(iB) + nobmS(iB));    % percent correct in go trials
    C.condNoGoPerfS(iB) = nobcS(iB) / (nobfS(iB) + nobcS(iB));    % percent correct in no-go trials
    C.condActionPerfS(iB) = nobhS(iB) / (nobhS(iB) + nobfS(iB));    % percent correct in trials when the mouse responded
    C.condRespS(iB) = (nobhS(iB) + nobfS(iB)) / nobS(iB);    % percent response in all trials
    C.condGoRespS(iB) = C.condGoPerfS(iB);    % percent response in go trials
    C.condNoGoRespS(iB) = 1 - C.condNoGoPerfS(iB);    % percent response in no-go trials
    C.condCorrectRespS(iB) = nobhS(iB) / (nobhS(iB) + nobcS(iB));    % percent response in correct trials
    C.condDiscriminationS(iB) = C.condGoRespS(iB) - C.condNoGoRespS(iB);   % hit - false alarm
    C.condDPrimeS(iB) = norminv(C.condGoRespS(iB)) - norminv(C.condNoGoRespS(iB));   % d' (SDT measure of discrimnability)
        
    nobL(iB) = length(intersect(find(bins),lowFR));   % number of trials in the bin
    binhits = intersect(find(ishit&bins),lowFR);
    nobhL(iB) = length(binhits);   % number of Hit trials in the bin
    binfas = intersect(find(isfa&bins),lowFR);
    nobfL(iB) = length(binfas);   % number of False alarm trials in the bin
    binmisss = intersect(find(ismiss&bins),lowFR);
    nobmL(iB) = length(binmisss);   % number of Miss trials in the bin
    bincrs = intersect(find(iscr&bins),lowFR);
    nobcL(iB) = length(bincrs);   % number of Correct rejection trials in the bin
    C.condPerfL(iB) = (nobhL(iB) + nobcL(iB)) / nobL(iB);    % percent correct
    C.condGoPerfL(iB) = nobhL(iB) / (nobhL(iB) + nobmL(iB));    % percent correct in go trials
    C.condNoGoPerfL(iB) = nobcL(iB) / (nobfL(iB) + nobcL(iB));    % percent correct in no-go trials
    C.condActionPerfL(iB) = nobhL(iB) / (nobhL(iB) + nobfL(iB));    % percent correct in trials when the mouse responded
    C.condRespL(iB) = (nobhL(iB) + nobfL(iB)) / nobL(iB);    % percent response in all trials
    C.condGoRespL(iB) = C.condGoPerfL(iB);    % percent response in go trials
    C.condNoGoRespL(iB) = 1 - C.condNoGoPerfL(iB);    % percent response in no-go trials
    C.condCorrectRespL(iB) = nobhL(iB) / (nobhL(iB) + nobcL(iB));    % percent response in correct trials
    C.condDiscriminationL(iB) = C.condGoRespL(iB) - C.condNoGoRespL(iB);   % hit - false alarm
    C.condDPrimeL(iB) = norminv(C.condGoRespL(iB)) - norminv(C.condNoGoRespL(iB));   % d' (SDT measure of discrimnability)
end