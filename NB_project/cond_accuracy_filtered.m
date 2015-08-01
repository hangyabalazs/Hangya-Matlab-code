function cond_accuracy_filtered
%COND_ACCURACY_FILTERED   Accuracy conditioned on reaction time.
%   COND_ACCURACYFILTERED calculates the percentage of hits relative to the
%   sum of hits and false alarms for specific reaction time intervals.
%   Sessions are prefiltered to remove 'impulsive' trials. Reaction times
%   are binned; hits and false alarms are partitioned according to the
%   reaction time bins. Sessions given in the input section are pooled
%   together. Conditional accuracy is plotted for different stimulus
%   intensities separately as well as pooled.
%
%   See also AUDITORY_GONOGO_PSYCHPLOT2 and IMPULSIVITY_FILTER.

%   Edit log: BH 7/8/11

% Input variables
animalID = 'n013';
sessionIDs = {'110629a','110630a','110701a','110704a','110705a'};
NUMsessions = length(sessionIDs);

% Preallocate variables
PerfGo = cell(1,NUMsessions);
PerfNoGo = cell(1,NUMsessions);
TimeGo = cell(1,NUMsessions);
TimeNoGo = cell(1,NUMsessions);
Sounds = cell(1,NUMsessions);
GoRTs = cell(1,NUMsessions);
NoGoRTs = cell(1,NUMsessions);
GoRTHit = cell(1,NUMsessions);
GoRTFA = cell(1,NUMsessions);

% Calculate performance for every session
for iS = 1:NUMsessions
    sessionID = sessionIDs{iS};
    
    % Load trial events structure
    fullpth = [getpref('cellbase','datapath') '\' animalID '\' sessionID '\'];
    TE = load([fullpth 'TrialEvents.mat']);
    TE = killnan(TE);   % remove NaNs
    TE0 = TE;
    [valid_trials, TE] = impulsivity_filter(TE);   % remove 'impulsive' trials
    
    % Calculate performance
    if isfield(TE0,'LightStimulation2') && ~all(isnan(TE0.LightStimulation2))
        lightoff_trial = find(TE0.LightStimulation2==0);
        [PerfGo{iS} PerfNoGo{iS} TimeGo{iS} TimeNoGo{iS}] = auditory_gonogo_psychplot2(animalID,sessionID,[],intersect(lightoff_trial,valid_trials));
    else
        [PerfGo{iS} PerfNoGo{iS} TimeGo{iS} TimeNoGo{iS}] = auditory_gonogo_psychplot2(animalID,sessionID,[],valid_trials);
    end
    
    % Find unique sound intensities
    [Sounds{iS}, ~, Ib] = unique(TE.StimulusDuration);
    NUMsounds = length(Sounds{iS});
    GoRTs{iS} = cell(1,NUMsounds);
    GoRTHit{iS} = cell(1,NUMsounds);
    GoRTFA{iS} = cell(1,NUMsounds);
    
    % Get all reaction times
    for iSound = 1:NUMsounds
        GoRTs{iS}{iSound} = TE.GoRT(Ib==iSound);   % RTs for a perticular sound intensity
        NoGoRTs{iS}{iSound} = TE.NoGoRT(Ib==iSound);
        GoRTHit{iS}{iSound} = nan2zero(TE.Hit(Ib==iSound));
        GoRTFA{iS}{iSound} = nan2zero(TE.FalseAlarm(Ib==iSound));
    end
    
end    % iS

% Get all sound intensities used in any of the sessions
diffsounds = unique(cell2mat(Sounds));
NUMdiffsounds = length(diffsounds);
GoRTs2 = cell(1,NUMdiffsounds);
NoGoRTs2 = cell(1,NUMdiffsounds);
GoRTHit2 = cell(1,NUMdiffsounds);
GoRTFA2 = cell(1,NUMdiffsounds);
GoRTHit3 = cell(1,NUMdiffsounds);
GoRTFA3 = cell(1,NUMdiffsounds);
condPerf = cell(1,NUMdiffsounds);

% Hit and miss conditioned on RT
for iDS = 1:NUMdiffsounds
    
    % Pool sessions for all different sound intensities
    for iS = 1:NUMsessions
        soundinx = find(Sounds{iS}==diffsounds(iDS));
        if ~isempty(soundinx)
            GoRTs2{iDS} = [GoRTs2{iDS} GoRTs{iS}{soundinx}];   % contains all RTs sorted according to different sounds - sessions are pooled
            NoGoRTs2{iDS} = [NoGoRTs2{iDS} NoGoRTs{iS}{soundinx}];
            GoRTHit2{iDS} = [GoRTHit2{iDS} GoRTHit{iS}{soundinx}];
            GoRTFA2{iDS} = [GoRTFA2{iDS} GoRTFA{iS}{soundinx}];
        end
    end
    
    % Histogram of RTs
    edges = (0:0.05:0.6);
    cnts = (edges(1:end-1) + edges(2:end)) / 2;
    [nmGoRTs binGoRTs] = histc(GoRTs2{iDS},edges);
    [nmNoGoRTs binNoGoRTs] = histc(NoGoRTs2{iDS},edges);
    
    % Partition performance according to different RT values
    NUMbins = length(cnts);
    GoRTHit3{iDS} = cell(1,NUMbins);
    GoRTFA3{iDS} = cell(1,NUMbins);
    condPerf{iDS} = zeros(1,NUMbins);
    for iB = 1:NUMbins
        goinx = binGoRTs == iB;
        GoRTHit3{iDS}{iB} = GoRTHit2{iDS}(goinx);   % partitioned according to RT values
        nogoinx = binNoGoRTs == iB;
        GoRTFA3{iDS}{iB} = GoRTFA2{iDS}(nogoinx);
        
        % Performance conditioned on RT
        condPerf{iDS}(iB) = sum(GoRTHit3{iDS}{iB}) / (sum(GoRTHit3{iDS}{iB}) + sum(GoRTFA3{iDS}{iB}));
    end
end

% Plot performance conditioned on RT
figure
for iDS = 1:NUMdiffsounds
    plot(cnts,condPerf{iDS},'Color',[iDS/NUMdiffsounds 0 0])
    hold on
end

% Regardless of sound intensities
condPerf2 = zeros(1,NUMbins);
GoRTHit4 = cell(1,NUMbins);
GoRTFA4 = cell(1,NUMbins);
for iB = 1:NUMbins
    for iDS = 1:NUMdiffsounds
        GoRTHit4{iB} = [GoRTHit4{iB} GoRTHit3{iDS}{iB}];   % pool different sounds together
        GoRTFA4{iB} = [GoRTFA4{iB} GoRTFA3{iDS}{iB}];
    end
    condPerf2(iB) = sum(GoRTHit4{iB}) / (sum(GoRTHit4{iB}) + sum(GoRTFA4{iB}));
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