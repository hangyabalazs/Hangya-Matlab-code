function fixTE2
%FIXTE2   Fix TrialEvent files.
%   FIXTE2 introduces a new event: 'FirstRestart' - the end of the first
%   ITI if there was a restart.
%
%   See also MAKETRIALEVENTS2_GONOGO

% Sessions
sessionids = listtag('session');

% Update 'TrialEvents'
NumSess = length(sessionids);   % number of sessions
cellbase_path = getpref('cellbase','datapath');   % CellBase location
for iC = 1:NumSess   % loop through all cells
    animalID = sessionids{iC,1};   % current animal
    sessionID = sessionids{iC,2};   % current session
    cellids = findcell('rat',animalID,'session',sessionID);
    cellid = cellids{1};   % choose one cell from the session to load TrialEvents
    if ~mod(iC,50)
        disp([num2str(iC) '   ' cellid])
    end
    try
        TE = loadcb(cellid,'TrialEvents');
    catch
        disp([cellid ': No TrialEvents file.'])
        continue
    end
    NumTrials = length(TE.TrialStart);   % number of trials
    [animal session] = cellid2tags(cellid);   % mouse ID and session ID
    sessionpath = fullfile(cellbase_path,animal,session);   % session full path   
        
    % First restart
    restartinx = cellfun(@(s)length(s)>1,TE.ITIBegins);   % more than one ITIs
    restarts = cellfun(@(s)s(1),TE.ITIEnds(restartinx));
    TE.FirstRestart = nan(1,NumTrials);
    TE.FirstRestart(restartinx) = restarts;
    
    % Save
    save(fullfile(sessionpath,'TrialEvents.mat'),'-struct','TE')
end