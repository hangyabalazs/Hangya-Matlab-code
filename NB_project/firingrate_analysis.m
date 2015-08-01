function FR = firingrate_analysis(cellid,varargin)
%FIRINGRATE_ANALYSIS   Firing rate.
%   FR = FIRINGRATE_ANALYSIS(CELLID) calculates firing rate in different
%   windows aligned to cue and feedback, in the full trial, before and
%   after the session is stopped.
%
%   See also REGRESSION_ANALYSIS.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   27-Aug-2012

%   Edit log: BH, 8/27/13

% Default arguments
prs = inputParser;
addRequired(prs,'cellid',@iscellid)
parse(prs,cellid,varargin{:});
g = prs.Results;

% Load trial events
try
    TE = loadcb(cellid,'TrialEvents');   % load events
    ST = loadcb(cellid,'EVENTSPIKES');   % load prealigned spikes
catch ME
    disp('There was no behavioral protocol for this session.')
    error(ME.message)
end

% Load stimulus events
try
    SE = loadcb(cellid,'StimEvents');   % load events
    SS = loadcb(cellid,'STIMSPIKES');   % load prealigned spikes
catch ME
    disp('There was no stimulus protocol for this session.')
%     error(ME.message)
end

% Load absoulte spike times
SP = loadcb(cellid);   % absoulte spike times

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
    fr_iti fr_iti05 fr_iti10 fr_trial] = deal(nan(NumTrials-1,1));
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
    
    % Firing rate in the full trial
    lim1 = TE.TrialStart(iT);
    lim2 = TE.TrialStart(iT+1);
    lspikes = SP;   % absolute spikes
    lspikes2 = lspikes(lspikes>lim1&lspikes<lim2);   % spikes in the time window
    fr_trial(iT) = length(lspikes2) / (lim2 - lim1);   % firing rate in the time window
end

% Firing rate before the session
lim1 = SP(1);   % first spike
lim2 = min(TE.TrialStart(1),SE.PulseOn(1));  % first trial or start of tagging
seg_presess = lim2 - lim1;  % segment length (in sec.)
presess = SP(SP>lim1&SP<lim2);   % post-last-hit spikes before session end
fr_presess = length(presess) / seg_presess;   % firing rate after the last hit

% Firing rate after stimulation, before session - tagging before behav.
lim1 = SE.PulseOn(find(SE.PulseOn<TE.TrialStart(1),1,'last')) + 0.1;  % stimulation session end
if ~isempty(lim1)
    lim2 = TE.TrialStart(find(TE.TrialStart>lim1,1,'first')); % first trial after stim.
    seg_poststim1 = lim2 - lim1;  % segment length (in sec.)
    poststim1 = SP(SP>lim1&SP<lim2);   % post-stimulation spikes before session end
    fr_poststim1 = length(poststim1) / seg_poststim1;   % firing rate after stimulation
else  % no stim. before the session
    seg_poststim1 = 0;
    fr_poststim1 = NaN;
end

% Firing rate before and after the last trial
lim1 = TE.TrialStart(find(~isnan(TE.TrialStart),1,'last'));  % last trial start (note that the last record can be NaN)
lim2 = SE.PulseOn(find(SE.PulseOn>lim1,1,'first'));  % stimulation session start
if isempty(lim2)
    lim2 = SP(end);   % if no stim.                                                        session after behav., use last spike
end
seg_postsess = lim2 - lim1;  % segment length (in sec.)
postsess = SP(SP>lim1&SP<lim2);   % post-session spikes
prepostsess = SP(SP>lim1-(lim2-lim1)&SP<lim1);   % spikes in a pre-session-end window of the same size
fr_postsess = length(postsess) / seg_postsess;  % firing rate after the end of the session
fr_prepostsess = length(prepostsess) / seg_postsess;   % firing rate before the end of the session

% Firing rate within the session after the last hit
lasthit = TE.TrialStart(find(TE.Hit==1,1,'last'));
lim1 = lasthit;   % last hit
lim2 = TE.TrialStart(find(~isnan(TE.TrialStart),1,'last'));  % last trial start (note that the last record can be NaN)
seg_stop = lim2 - lim1;  % segment length (in sec.)
stopspikes = SP(SP>lim1&SP<lim2);   % post-last-hit spikes before session end
fr_stop = length(stopspikes) / (lim2 - lim1);   % firing rate after the last hit

% Firing rate after last stimulation - tagging after behaviour
lim1 = SE.PulseOn(end) + 0.1;  % stimulation session end
lim2 = SP(end);   % last spike
seg_poststim2 = lim2 - lim1;  % segment length (in sec.)
poststim2 = SP(SP>lim1&SP<lim2);   % post-stimulation spikes before session end
fr_poststim2 = length(poststim2) / seg_poststim2;   % firing rate after stimulation

% Position data
[animalID sessionID] = cellid2tags(cellid);
fullpth = fullfile(getpref('cellbase','datapath'),animalID,sessionID);
fn = fullfile(fullpth,'VT1.nvt');
[TimeStamps, ExtractedX, ExtractedY, ExtractedAngle, Targets, Points, NlxHeader] =...
    Nlx2MatVT(fn,[1 1 1 1 1 1],1,1);
TimeStamps = TimeStamps / 1e6;

% Load LFP
cscname = fullfile(fullpth,'CSC5.ncs');   % filename for the Neuralynx CSC file
[TimeStamps_CSC, ChanNum, SampleFrequency, NumValSamples, Samples, NlxHeader] = ...
    Nlx2MatCSC(cscname,[1 1 1 1 1],1,1,1);  % mPFC LFP
TimeStamps_CSC = TimeStamps_CSC / 1e6;
Samples = Samples(1:32:end);
LFP = -1 * Samples(:);
sr = 16 / (mean(diff(TimeStamps_CSC)));
dt = 1 /sr;
time = TimeStamps_CSC(1):dt:TimeStamps_CSC(1)+dt*(length(LFP)-1);
time_orig = repmat(TimeStamps_CSC,16,1) + repmat((0:15)'*dt,1,length(TimeStamps_CSC));
time_orig = time_orig(:)';

% Output
FR.stim_0_05 = fr_stim_0_05;
FR.stim_0_10 = fr_stim_0_10;
FR.stim2resp = fr_stim2resp;
FR.resp_0_05 = fr_resp_0_05;
FR.resp_0_10 = fr_resp_0_10;
FR.resp_02_0 = fr_resp_02_0;
FR.stim_025_0 = fr_stim_025_0;
FR.stim_05_0 = fr_stim_05_0;
FR.stim_05_02 = fr_stim_05_02;
FR.stim_10_0 = fr_stim_10_0;
FR.iti = fr_iti;
FR.iti05 = fr_iti05;
FR.iti10 = fr_iti10;
FR.trial = fr_trial;
FR.presess = fr_presess;
FR.prepostsess = fr_prepostsess;
FR.postsess = fr_postsess;
FR.stop = fr_stop;
FR.poststim1 = fr_poststim1;
FR.poststim2 = fr_poststim2;

% keyboard

figure;plot(fr_trial)
fr_trial2 = fr_trial;
fr_trial2(TE.Hit==1|TE.FalseAlarm==1) = NaN;
figure;plot(fr_trial2,'.')
figure;plot(fr_stim_05_0)
nanmean(fr_trial)
nanmean(fr_trial(TE.CorrectRejection==1))
cellid
% !start AcroRd32.exe g:\NB_cellbase\n023\111220a\n023_111220a_clusters.pdf
figure;plot(farate)
figure;plot(hitrate)
figure;plot(ExtractedX,ExtractedY,'.')
figure;plot(ExtractedX(TimeStamps>TE.TrialStart(end-1)),ExtractedY(TimeStamps>TE.TrialStart(end-1)),'.')
figure;plot(ExtractedX(TimeStamps>SE.PulseOn(1)&TimeStamps<SE.PulseOn(end)),ExtractedY(TimeStamps>SE.PulseOn(1)&TimeStamps<SE.PulseOn(end)),'.')
figure;plot(ExtractedX(TimeStamps<TE.TrialStart(1)),ExtractedY(TimeStamps<TE.TrialStart(1)),'.')
% figure;plot(LFP)

keyboard