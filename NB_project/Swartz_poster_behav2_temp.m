function Swartz_poster_behav2

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
GoRTHit = cell(1,NUMsessions);
GoRTMiss = cell(1,NUMsessions);

% Calculate performance for every session
for iS = 1:NUMsessions
    sessionID = sessionIDs{iS};
    
    % Load trial events structure
    fullpth = [getpref('cellbase','datapath') '\' animalID '\' sessionID '\'];
    TE = load([fullpth 'TrialEvents.mat']);
    TE = killnan(TE);   % remove NaNs
    
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
    GoRTs{iS} = cell(1,NUMsounds);
    GoRTHit{iS} = cell(1,NUMsounds);
    GoRTMiss{iS} = cell(1,NUMsounds);
    
    % Get all reaction times
    for iSound = 1:NUMsounds
        GoRTs{iS}{iSound} = TE.GoRT(Ib==iSound);   % RTs for a perticular sound intensity
        GoRTHit{iS}{iSound} = nan2zero(TE.Hit(Ib==iSound));
        GoRTMiss{iS}{iSound} = nan2zero(TE.Miss(Ib==iSound));
    end
    
end    % iS

% Get all sound intensities used in any of the sessions
diffsounds = unique(cell2mat(Sounds));
NUMdiffsounds = length(diffsounds);
GoRTs2 = cell(1,NUMdiffsounds);
GoRTHit2 = cell(1,NUMdiffsounds);
GoRTMiss2 = cell(1,NUMdiffsounds);
GoRTHit3 = cell(1,NUMdiffsounds);
GoRTMiss3 = cell(1,NUMdiffsounds);
condPerf = cell(1,NUMdiffsounds);

% Hit and miss conditioned on RT
for iDS = 1:NUMdiffsounds
    
    % Pool sessions for all different sound intensities
    for iS = 1:NUMsessions
        soundinx = find(Sounds{iS}==diffsounds(iDS));
        if ~isempty(soundinx)
            GoRTs2{iDS} = [GoRTs2{iDS} GoRTs{iS}{soundinx}];   % contains all RTs sorted according to different sounds - sessions are pooled
            GoRTHit2{iDS} = [GoRTHit2{iDS} GoRTHit{iS}{soundinx}];
            GoRTMiss2{iDS} = [GoRTMiss2{iDS} GoRTMiss{iS}{soundinx}];
        end
    end
    
    % Histogram of RTs
    edges = (0.1:0.05:0.4);
    cnts = (edges(1:end-1) + edges(2:end)) / 2;
    [nmGoRTs binGoRTs] = histc(GoRTs2{iDS},edges);
    
    % Partition performance according to different RT values
    NUMbins = length(cnts);
    GoRTHit3{iDS} = cell(1,NUMbins);
    GoRTMiss3{iDS} = cell(1,NUMbins);
    condPerf{iDS} = zeros(1,NUMbins);
    for iB = 1:NUMbins
        inx = binGoRTs == iB;
        GoRTHit3{iDS}{iB} = GoRTHit2{iDS}(inx);   % partitioned according to RT values
        GoRTMiss3{iDS}{iB} = GoRTMiss2{iDS}(inx);
        
        % Performance conditioned on RT
        condPerf{iDS}(iB) = sum(TE.GoRTHit3{iDS}{iB}) / (sum(GoRTHit3{iDS}{iB}) + sum(GoRTMiss3{iDS}{iB}));
    end
end

% Plot performance conditioned on RT
figure
for iDS = 1:NUMdiffsounds
    plot(cnts,condPerf{iDS},'Color',[iDS/NUMdiffsounds iDS/NUMdiffsounds iDS/NUMdiffsounds])
    hold on
end

% -------------------------------------------------------------------------
function TE = killnan(TE)

% Find NaNs
naninx = isnan(TE.SoundDuration);

% Remove NaNs
fnm = fieldnames(TE);
for k = 1:length(fieldnames(TE))
    TE.(fnm{k})(naninx) = [];
end