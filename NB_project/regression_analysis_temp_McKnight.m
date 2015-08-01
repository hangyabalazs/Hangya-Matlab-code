function [R p] = regression_analysis_temp_McKnight(cellid,varargin)
%REGRSSION_ANALYSIS   Linear regression for firing rate.
%   REGRESSION_ANALYSIS(CELLID) calculates linear regression coefficients
%   for firing rate variables (aligned to stimulus and response onset) and
%   various task-related regressors.
%
%   See also NBREGRESSION.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   6-Oct-2012

%   Edit log: 10/6/12, 7/30/13

% ideas: number of restart ITIs; blockwise 'difficulty'; difficulty in
% previous trial

% Default arguments
prs = inputParser;
addRequired(prs,'cellid',@iscellid)
g = parse(prs,cellid,varargin{:});

% Load trial events
try
    TE = loadcb(cellid,'TrialEvents');   % load events
    ST = loadcb(cellid,'EVENTSPIKES');   % load prealigned spikes
catch ME
    disp('There was no behavioral protocol for ths session.')
    error(ME.message)
end

% Checking whether 'DeliverFeedback' event is available
sesstype = getvalue('session_type',cellid);
if isequal(sesstype,{'feedbackdelay'})
    alignevent_fa = 'DeliverFeedback';
    alignevent_hit = 'DeliverFeedback';
else
    alignevent_fa = 'LeftPortIn';
    alignevent_hit = 'LeftWaterValveOn';
end

% Relative spike times
stim_pos = findcellstr(ST.events(:,1),'StimulusOn');   % tone onset
response_pos_fa = findcellstr(ST.events(:,1),alignevent_fa);   % animal's respones, false alarms
response_pos_hit = findcellstr(ST.events(:,1),alignevent_hit);   % animal's respones, hits 
spikes_stimulus = ST.event_stimes{stim_pos};   % spikes relative to tone onset
spikes_response = cell(size(spikes_stimulus));
spikes_response(TE.FalseAlarm==1) = ST.event_stimes{response_pos_fa}(TE.FalseAlarm==1);   % spikes relative to response onset, false alarms
spikes_response(TE.Hit==1) = ST.event_stimes{response_pos_hit}(TE.Hit==1);   % spikes relative to response onset, hits

% Calculate regressors and dependent variables
NumTrials = length(TE.TrialStart);   % number of trials
wn = 20;   % window size (number of trials)
[hits,fas,gos,nogos,hitrate,farate,discrimination,dprime,engagement,...
    correct,incorrect,responded,skipped,...
    gort,nogort,gortint,nogortint,gortloc,nogortloc,loudness,iti,...
    fr_stim_0_05,fr_stim_0_10,fr_stim2resp,fr_resp_0_05,fr_resp_0_10,fr_resp_02_0,...
    fr_stim_025_0, fr_stim_05_0, fr_stim_05_02, fr_stim_10_0,...
    fr_iti fr_iti05 fr_iti10] = deal(nan(NumTrials-1,1));
for iT = 1:NumTrials-1     % last trial may be incomplete
    
    % Outcome
    if iT >= wn
        hits(iT) = sum(~isnan(TE.Hit(iT-wn+1:iT)));   % number of hits
        fas(iT) = sum(~isnan(TE.FalseAlarm(iT-wn+1:iT)));   % number of false alarms
        gos(iT) = sum(~isnan(TE.Hit(iT-wn+1:iT))|~isnan(TE.Miss(iT-wn+1:iT)));   % number of go tones
        nogos(iT) = sum(~isnan(TE.FalseAlarm(iT-wn+1:iT))|~isnan(TE.CorrectRejection(iT-wn+1:iT)));   % number of no-go tones
    end
    
    hitrate(iT) = hits(iT) / gos(iT);   % hit rate
    farate(iT) = fas(iT) / nogos(iT);   % false alarm rate
    discrimination(iT) = hitrate(iT) - farate(iT);   % hit - false alarm
    dprime(iT) = norminv(hitrate(iT)) - norminv(farate(iT));   % d' (SDT measure of discrimnability)
    engagement(iT) = hitrate(iT) + farate(iT);   % hit + false alarm
    
    correct(iT) = isequal(TE.Hit(iT),1) | isequal(TE.CorrectRejection(iT),1);
    incorrect(iT) = isequal(TE.FalseAlarm(iT),1) | isequal(TE.Miss(iT),1);
    responded(iT) = isequal(TE.Hit(iT),1) | isequal(TE.FalseAlarm(iT),1);
    skipped(iT) = isequal(TE.Miss(iT),1) | isequal(TE.CorrectRejection(iT),1);
    
    % Reaction time
    gort(iT) = TE.GoRT(iT);   % reaction time
    nogort(iT) = TE.NoGoRT(iT);   % response time for false alarms
    if iT >= wn
        gortint(iT) = nanmean(TE.GoRT(iT-wn+1:iT));   % reaction time averaged over the window
        nogortint(iT) = nanmean(TE.NoGoRT(iT-wn+1:iT));   % no-go response time averaged over the window
    end
    gortloc(iT) = gort(iT) - gortint(iT);   % current reaction time relative to window average
    nogortloc(iT) = nogort(iT) - nogortint(iT);   % current no-go response time relative to window average 
    
    % Sound intensity
    loudness(iT) = TE.StimulusDuration(iT);   % tone intensity
    
    % Foreperiod (expectancy)
    iti(iT) = TE.ITIDistribution(iT);   % length of foreperiod (s)
        
    % Firing rate 0.5 s after stimulus onset
    lim1 = 0;  % start of time window
    lim2 = 0.5;   % end of time window
    lspikes = spikes_stimulus{iT};   % spikes relative to the event
    lspikes2 = lspikes(lspikes>lim1&lspikes<lim2);   % spikes in the time window
    fr_stim_0_05(iT) = length(lspikes2) / (lim2 - lim1);   % firing rate in the time window
    
    % Firing rate 1 s after stimulus onset
    lim1 = 0;  % start of time window
    lim2 = 1;   % end of time window
    lspikes = spikes_stimulus{iT};   % spikes relative to the event
    lspikes2 = lspikes(lspikes>lim1&lspikes<lim2);   % spikes in the time window
    fr_stim_0_10(iT) = length(lspikes2) / (lim2 - lim1);   % firing rate in the time window
    
    % Firing rate in from stimulus to response
    lim1 = 0;
    lim2 = TE.ReactionTime(iT);
    lspikes = spikes_stimulus{iT};   % spikes relative to the event
    lspikes2 = lspikes(lspikes>lim1&lspikes<lim2);   % spikes in the time window
    fr_stim2resp(iT) = length(lspikes2) / (lim2 - lim1);   % firing rate in the time window
    
    % Firing rate 0.25 s before stimulus onset
    lim1 = -0.25;  % start of time window
    lim2 = 0;   % end of time window
    lspikes = spikes_stimulus{iT};   % spikes relative to the event
    lspikes2 = lspikes(lspikes>lim1&lspikes<lim2);   % spikes in the time window
    fr_stim_025_0(iT) = length(lspikes2) / (lim2 - lim1);   % firing rate in the time window
    
    % Firing rate 0.5 s before stimulus onset
    lim1 = -0.5;  % start of time window
    lim2 = 0;   % end of time window
    lspikes = spikes_stimulus{iT};   % spikes relative to the event
    lspikes2 = lspikes(lspikes>lim1&lspikes<lim2);   % spikes in the time window
    fr_stim_05_0(iT) = length(lspikes2) / (lim2 - lim1);   % firing rate in the time window
    
    % Firing rate 0.5-0.2 s before stimulus onset
    lim1 = -0.5;  % start of time window
    lim2 = -0.2;   % end of time window
    lspikes = spikes_stimulus{iT};   % spikes relative to the event
    lspikes2 = lspikes(lspikes>lim1&lspikes<lim2);   % spikes in the time window
    fr_stim_05_02(iT) = length(lspikes2) / (lim2 - lim1);   % firing rate in the time window
    
    % Firing rate 1 s before stimulus onset
    lim1 = -1;  % start of time window
    lim2 = 0;   % end of time window
    lspikes = spikes_stimulus{iT};   % spikes relative to the event
    lspikes2 = lspikes(lspikes>lim1&lspikes<lim2);   % spikes in the time window
    fr_stim_10_0(iT) = length(lspikes2) / (lim2 - lim1);   % firing rate in the time window
    
    % Firing rate 0.5 s after response onset
    lim1 = 0;  % start of time window
    lim2 = 0.5;   % end of time window
    lspikes = spikes_response{iT};   % spikes relative to the event
    lspikes2 = lspikes(lspikes>lim1&lspikes<lim2);   % spikes in the time window
    fr_resp_0_05(iT) = length(lspikes2) / (lim2 - lim1);   % firing rate in the time window
    
    % Firing rate 1 s after response onset
    lim1 = 0;  % start of time window
    lim2 = 1;   % end of time window
    lspikes = spikes_response{iT};   % spikes relative to the event
    lspikes2 = lspikes(lspikes>lim1&lspikes<lim2);   % spikes in the time window
    fr_resp_0_10(iT) = length(lspikes2) / (lim2 - lim1);   % firing rate in the time window
    
    % Firing rate 0.2 s before response onset
    lim1 = -0.2;  % start of time window
    lim2 = 0;   % end of time window
    lspikes = spikes_response{iT};   % spikes relative to the event
    lspikes2 = lspikes(lspikes>lim1&lspikes<lim2);   % spikes in the time window
    fr_resp_02_0(iT) = length(lspikes2) / (lim2 - lim1);   % firing rate in the time window
    
    % Firing rate in the foreperiod
    itiwin = TE.ITIDistribution(iT);   % ITI
    lim1 = -itiwin;
    lim2 = 0;
    lspikes = spikes_stimulus{iT};   % spikes relative to the event
    lspikes2 = lspikes(lspikes>lim1&lspikes<lim2);   % spikes in the time window
    fr_iti(iT) = length(lspikes2) / (lim2 - lim1);   % firing rate in the time window
    
    % Firing rate in the foreperiod, but max 0.5 s
    itiwin = TE.ITIDistribution(iT);   % ITI
    lim1 = max(-itiwin,-0.5);
    lim2 = 0;
    lspikes = spikes_stimulus{iT};   % spikes relative to the event
    lspikes2 = lspikes(lspikes>lim1&lspikes<lim2);   % spikes in the time window
    fr_iti05(iT) = length(lspikes2) / (lim2 - lim1);   % firing rate in the time window
    
    % Firing rate in the foreperiod, but max 1 s
    itiwin = TE.ITIDistribution(iT);   % ITI
    lim1 = max(-itiwin,-1);
    lim2 = 0;
    lspikes = spikes_stimulus{iT};   % spikes relative to the event
    lspikes2 = lspikes(lspikes>lim1&lspikes<lim2);   % spikes in the time window
    fr_iti10(iT) = length(lspikes2) / (lim2 - lim1);   % firing rate in the time window
end

% Regression
hitinx = find(~isnan(TE.FalseAlarm));
% hitinx = intersect(find(~isnan(TE.Hit)),(wn:NumTrials-1));
% hitinx = intersect(intersect(find(~isnan(TE.Hit)),(wn:NumTrials-1)),find(~isnan(dprime)&~ismember(dprime,[Inf,-Inf])));
% hitinx = intersect(hitinx,find(gort>0.15));
y = loudness(hitinx);
X = [ones(length(hitinx),1) fr_stim_0_05(hitinx)];
[b,bint,r,rint,stats] = regress(y,X);
p.stim2resp = stats(3);
pR = corrcoef(y,X(:,2));
R.stim2resp = pR(3);