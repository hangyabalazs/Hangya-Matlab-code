function session_type = nbsessionclass
%NBSESSIONCLASS   Classify behavioral sessions.
%   SESSION_TYPE = NBSESSIONCLASS classifies sessions based on the behavior
%   file to 'no behavior' (no behavior file), 'feedbackdelay' (variable
%   delay between response and feedback), 'pulsepal' (laser stimulation)
%   and 'gonogo' (basic go-nogo task) sessions (SESSION_TYPE).

% Load CellBase
loadcb

% Session tags
NumCells = length(CELLIDLIST);   % number of cells in CellBase
session_type = cell(NumCells,1);
for iC = 1:NumCells   % loop through all cells
    cellid = CELLIDLIST(iC);   % current cell ID
    fullpth = [cellid2fnames(cellid,'session') filesep];   % pathname for the session
    
    % Find behavior file
    fls = dir(fullpth);
    fls = {fls.name};  % filenames in deirectory
    inx = find(strncmp(fls,'data_@auditory',14));
        
    % Session type
    if isempty(inx)
        session_type{iC} = 'no behavior';   % no behav. file was found
        continue
    else
        behav_file = fls{inx};  % behavior file
    end
    if ~isempty(strfind(behav_file,'feedbackdelay'))
        session_type{iC} = 'feedbackdelay';   % session w feedback-delay
    elseif ~isempty(strfind(behav_file,'pulsepal'))
        session_type{iC} = 'pulsepal';   % session w PulsePal stimulation
    elseif ~isempty(strfind(behav_file,'pavlovian'))
        session_type{iC} = 'pavlovian';   % Pavlovian session
    else
        session_type{iC} = 'gonogo';   % basic auditory go-nogo
    end
end