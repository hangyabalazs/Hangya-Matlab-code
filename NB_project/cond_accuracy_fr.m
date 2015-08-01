function cond_accuracy_fr
%COND_ACCURACY_FR   Accuracy conditioned on firing rate.
%   COND_ACCURACY_FR calculates the percentage of hits relative to the sum
%   of hits and false alarms for specific firing rate intervals of a given
%   cell. Firing rates are binned; hits and false alarms are partitioned
%   according to the firing rate bins. Sessions given in the input
%   section are pooled together. Conditional accuracy is plotted for
%   different stimulus intensities separately as well as pooled.
%
%   See also COND_ACCURACY.

%   Edit log: BH 7/7/11

% Input variables
animalID = 'n013';
sessionIDs = {'110628a'};
NUMsessions = length(sessionIDs);
sessionID = sessionIDs{1};
fullpth = [getpref('cellbase','datapath') '\' animalID '\' sessionID '\'];
TE = load([fullpth 'TrialEvents.mat']);
load([fullpth 'EVENTSPIKES5_4.mat'])
spikes_stimon = event_stimes{1};

% Prestimulus frequency
NUMtrials = length(spikes_stimon);
prestimfreq = nan(1,NUMtrials);
for k = 1:NUMtrials
    lspikes = spikes_stimon{k};
    lspikes2 = lspikes(lspikes>-2&lspikes<0);   % time window: one sec before stimulus onset
    prestimfreq(k) = length(lspikes2) / 2;
end

% Preallocate variables
PerfGo = cell(1,NUMsessions);
PerfNoGo = cell(1,NUMsessions);
TimeGo = cell(1,NUMsessions);
TimeNoGo = cell(1,NUMsessions);
Sounds = cell(1,NUMsessions);
GoFRs = cell(1,NUMsessions);
NoGoFRs = cell(1,NUMsessions);
GoFRHit = cell(1,NUMsessions);
GoFRFA = cell(1,NUMsessions);

% Calculate performance for every session
for iS = 1:NUMsessions
    sessionID = sessionIDs{iS};
    
    % Load trial events structure
    fullpth = [getpref('cellbase','datapath') '\' animalID '\' sessionID '\'];
    TE = load([fullpth 'TrialEvents.mat']);
    TE = killnan(TE);   % remove NaNs
    
    % Response types
    ishit = logical(nan2zero(TE.Hit));
    ismiss = logical(nan2zero(TE.Miss));
    isfa = logical(nan2zero(TE.FalseAlarm));
    iscr = logical(nan2zero(TE.CorrectRejection));
    prestimfreq_hm = prestimfreq;
    prestimfreq_hm(~(ishit|ismiss)) = NaN;
%     prestimfreq_hm(TE.GoRT<0.1) = NaN;   % exclude 'impulsive' trials
%     prestimfreq_hm(TE.NoGoRT<0.1) = NaN;
    prestimfreq_fc = prestimfreq;
    prestimfreq_fc(~(isfa|iscr)) = NaN;
%     prestimfreq_fc(TE.GoRT<0.1) = NaN;   % exclude 'impulsive' trials
%     prestimfreq_fc(TE.NoGoRT<0.1) = NaN;
    
    % Calculate performance
    if isfield(TE,'LightStimulation2')
        lightoff_trial = find(TE.LightStimulation2==0);
        [PerfGo{iS} PerfNoGo{iS} TimeGo{iS} TimeNoGo{iS}] = auditory_gonogo_psychplot2(animalID,sessionID,[],lightoff_trial);
    else
        [PerfGo{iS} PerfNoGo{iS} TimeGo{iS} TimeNoGo{iS}] = auditory_gonogo_psychplot2(animalID,sessionID);
    end
    
    % Find unique sound intensities
    [Sounds{iS}, ~, Ib] = unique(TE.StimulusDuration);
    NUMsounds = length(Sounds{iS});
    GoFRs{iS} = cell(1,NUMsounds);
    GoFRHit{iS} = cell(1,NUMsounds);
    GoFRFA{iS} = cell(1,NUMsounds);
    
    % Get all reaction times
    for iSound = 1:NUMsounds
        GoFRs{iS}{iSound} = prestimfreq_hm(Ib==iSound);   % FRs for a perticular sound intensity
        NoGoFRs{iS}{iSound} = prestimfreq_fc(Ib==iSound);   % 'impulsive' (faster than 100 ms) responses excluded
        GoFRHit{iS}{iSound} = nan2zero(TE.Hit(Ib==iSound));
        GoFRFA{iS}{iSound} = nan2zero(TE.FalseAlarm(Ib==iSound));
    end
    
end    % iS

% Get all sound intensities used in any of the sessions
diffsounds = unique(cell2mat(Sounds));
NUMdiffsounds = length(diffsounds);
GoFRs2 = cell(1,NUMdiffsounds);
NoGoFRs2 = cell(1,NUMdiffsounds);
GoFRHit2 = cell(1,NUMdiffsounds);
GoFRFA2 = cell(1,NUMdiffsounds);
GoFRHit3 = cell(1,NUMdiffsounds);
GoFRFA3 = cell(1,NUMdiffsounds);
condPerf = cell(1,NUMdiffsounds);

% Hit and miss conditioned on FR
for iDS = 1:NUMdiffsounds
    
    % Pool sessions for all different sound intensities
    for iS = 1:NUMsessions
        soundinx = find(Sounds{iS}==diffsounds(iDS));
        if ~isempty(soundinx)
            GoFRs2{iDS} = [GoFRs2{iDS} GoFRs{iS}{soundinx}];   % contains all FRs sorted according to different sounds - sessions are pooled
            NoGoFRs2{iDS} = [NoGoFRs2{iDS} NoGoFRs{iS}{soundinx}];
            GoFRHit2{iDS} = [GoFRHit2{iDS} GoFRHit{iS}{soundinx}];
            GoFRFA2{iDS} = [GoFRFA2{iDS} GoFRFA{iS}{soundinx}];
        end
    end
    
    % Histogram of FRs
    edges = (10:10:40);
    cnts = (edges(1:end-1) + edges(2:end)) / 2;
    [nmGoFRs binGoFRs] = histc(GoFRs2{iDS},edges);
    [nmNoGoFRs binNoGoFRs] = histc(NoGoFRs2{iDS},edges);
    
    % Partition performance according to different FR values
    NUMbins = length(cnts);
    GoFRHit3{iDS} = cell(1,NUMbins);
    GoFRFA3{iDS} = cell(1,NUMbins);
    condPerf{iDS} = zeros(1,NUMbins);
    for iB = 1:NUMbins
        goinx = binGoFRs == iB;
        GoFRHit3{iDS}{iB} = GoFRHit2{iDS}(goinx);   % partitioned according to FR values
        nogoinx = binNoGoFRs == iB;
        GoFRFA3{iDS}{iB} = GoFRFA2{iDS}(nogoinx);
        
        % Performance conditioned on FR
        condPerf{iDS}(iB) = sum(GoFRHit3{iDS}{iB}) / (sum(GoFRHit3{iDS}{iB}) + sum(GoFRFA3{iDS}{iB}));
    end
end

% Plot performance conditioned on FR
figure
for iDS = 1:NUMdiffsounds
    plot(cnts,condPerf{iDS},'Color',[iDS/NUMdiffsounds 0 0])
    hold on
end

% Regardless of sound intensities
condPerf2 = zeros(1,NUMbins);
GoFRHit4 = cell(1,NUMbins);
GoFRFA4 = cell(1,NUMbins);
for iB = 1:NUMbins
    for iDS = 1:NUMdiffsounds
        GoFRHit4{iB} = [GoFRHit4{iB} GoFRHit3{iDS}{iB}];   % pool different sounds together
        GoFRFA4{iB} = [GoFRFA4{iB} GoFRFA3{iDS}{iB}];
    end
    condPerf2(iB) = sum(GoFRHit4{iB}) / (sum(GoFRHit4{iB}) + sum(GoFRFA4{iB}));
end

figure   % plot
plot(cnts,condPerf2)
keyboard

% -------------------------------------------------------------------------
function TE = killnan(TE)

% Find NaNs
naninx = isnan(TE.SoundDuration);

% Remove NaNs
fnm = fieldnames(TE);
for k = 1:length(fieldnames(TE))
    TE.(fnm{k})(naninx) = [];
end