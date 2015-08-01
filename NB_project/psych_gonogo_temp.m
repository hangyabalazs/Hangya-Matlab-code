function psych_gonogo(cellids)

% Input arguments
mode = 'restrict';   % include only those sessions that contain the input cell IDs

% Animals
mice = listtag('animal');
NumMice = length(mice);   % number of animals

% Performance
for iM = 1:NumMice   % loop through mice
    animalID = mice(iM)   % current mouse
    sessions = findsession('animal',animalID);   % sessions of the current mouse
    NumSessions = size(sessions,1);   % number of sessions
    
    % Per animal
    for iS = 1:NumSessions   % loop through sessions
        sessionID = sessions(iS)   % current session
        cellIDs = findcell('rat',animalID,'session',sessionID);   % cells of the current session
        if isequal(mode,'restrict') && isempty(intersect(cellids,cellIDs))
            continue   % skip session if it does not contain any of the input cell IDs (e.g. no valid cluster from NB)
        end
        
        % Load trial events
        cellid = cellids{1};   % choose one cell from the session to load TrialEvents
        try
            TE = loadcb(cellid,'TrialEvents');
        catch
            disp([cellid ': No TrialEvents file.'])
            continue
        end
        
        % Session performance
        if isfield(TE,'LightStimulation2')   % exclude trials with light-stimulation (affects only one session of NB CellBase)
            lightoff_trial = find(TE.LightStimulation2==0);
        else
            lightoff_trial = 1:length(TE.TrialStart);
        end
        [PerfGo2 PerfNoGo2 TimeGo2 TimeNoGo2] = ...
            auditory_gonogo_psychplot3(animalID,sessionID,[],lightoff_trial);
        
    end
    
end