function learning_curve
%COND_ACCURACY   Accuracy conditioned on reaction time.
%   COND_ACCURACY calculates the percentage of hits relative to the sum of
%   hits and false alarms for specific reaction time intervals. Reaction
%   times are binned; hits and false alarms are partitioned according to
%   the reaction time bins. Sessions given in the input section are pooled
%   together. Conditional accuracy is plotted for different stimulus
%   intensities separately as well as pooled.
%
%   See also AUDITORY_GONOGO_PSYCHPLOT2.

%   Edit log: BH 7/6/11

% Input variables
animalID = 'n013';
animalID2 = 'nb013';
inpdir = 'c:\Balazs\_data\NB\behav\nb013\';
psi = dir(inpdir);
sessionIDs = {psi(3:end).name};
NUMsessions = length(sessionIDs);

% Preallocate variables
PerfGo = zeros(1,NUMsessions);
PerfNoGo = zeros(1,NUMsessions);

% Calculate performance for every session
for iS = 1:NUMsessions
    sessionID = sessionIDs{iS}(end-10:end-4);
    
    % Load trial events structure
    fullpth = inpdir;
    TE = solo2trialevents2_auditory_gonogo([fullpth 'data_@auditory_gonogo_balazs_' animalID2 '_' sessionID '.mat'],0);
    TE = killnan(TE);   % remove NaNs; note that this TE is not synchronized to Neuralynx!
    
    % Performance
    PerfGo(iS) = nansum(TE.Hit) / (nansum(TE.Hit) + nansum(TE.Miss));
    PerfNoGo(iS) = nansum(TE.FalseAlarm) / (nansum(TE.CorrectRejection) + nansum(TE.FalseAlarm));
end

% Plot performance for all sessions
figure
plot(PerfGo,'g')
hold on
plot(PerfNoGo,'r')
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