function fixTE(cellids)
%FIXTE   Fix TrialEvent files.
%   FIXTE(CELLIDS) updates TrialEvent files for a set of cells (CELLIDS).
%   Some events were erroneously saved in absolute time stamps as opposed
%   to time stamps relative to 'TrialStart'. These time stamps became
%   meaningless after synchronization (see MAKETRIALEVENTS2_GONOGO). FIXTE
%   corrects these time stamps based on the original, un-synchronized
%   TrialEvents file. The affected fields are 'TrialEnd' (for all cells),
%   'ITIBegins', 'ITIEnds', NoPokeITIBegins', 'NoPokeITIEnds' (for initial
%   cells saved by an older version of 'solo2trialevents').
%
%   FIXTE also introduces two new events: 'LastITIBegins' and 'LastITIEnds'
%   - boundaries of the ITI after the last restart.
%
%   See also MAKETRIALEVENTS2_GONOGO

% Cell IDs
if nargin < 1 || isempty(cellids)
    loadcb
    cellids = CELLIDLIST;
else
    if isnumeric(cellids)
        loadcb
        cellids = CELLIDLIST(cellids);   % index set to CELLIDLIST
    elseif ischar(cellids)
        cellids = {cellids};   % only one cellID passed
    elseif iscellstr(cellids)
        % list of cell IDs
    else
        error('fixTE:inputArg','Unsupported format for cell IDs.')
    end
end

% Update 'TrialEvents'
NumCells = length(cellids);   % number of cells
cellbase_path = getpref('cellbase','datapath');   % CellBase location
for iC = 1:NumCells   % loop through all cells
    cellid = cellids{iC};   % current cell ID
    if ~mod(iC,100)
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
    fn = fullfile(sessionpath,'TE.mat');
    oldTE = load(fn);   % original (un-synchronized) TrialEvents file
    
    % Match TrialEvent structures
    so = TE.StimulusOn;   % StimulusOn should be the same
    oldso = oldTE.StimulusOn;
    binx = [];
    next = 0;
    while ~isequalwithequalnans(so,oldso)   % while there is a mismatch
        tinx = find(abs(so-oldso(1:length(so)))>0,1,'first');   % find first mismatch
        oldso(tinx) = [];
        binx = [binx tinx+next];   % collect bad indices
        next = next + 1;
    end
    if ~isempty(binx)
        ginx = setdiff(1:length(oldTE.TrialStart),binx);   % good indeices
        oldTE = shortenTE(oldTE,ginx);   % eliminate trials from the un-sync TrialEvents
    end
    if ~isequalwithequalnans(TE.StimulusOn,oldTE.StimulusOn)   % safety check
        error('fixTE:mismatch','Failed to match trials.')
    end
    
    % TrialEnd
    df = diff(TE.TrialEnd);
    df = df(~isnan(df));
    if all(df>0)   % monotonically increasing time stamps indicates that they are in abs. times
        TE.TrialEnd = TE.TrialEnd - oldTE.TrialStart;   % make it rel. to TrialStart
    end
    
    % ITIBegins
    df = diff(cellfun(@(s)s(1),TE.ITIBegins(cellfun(@(s)~isempty(s),TE.ITIBegins))));
    df = df(~isnan(df));
    if all(df>0)   % monotonically increasing time stamps indicates that they are in abs. times
        for iT = 1:NumTrials
            TE.ITIBegins{iT} = TE.ITIBegins{iT} - oldTE.TrialStart(iT);   % make it rel. to TrialStart
        end
    end
    
    % ITIEnds
    df = diff(cellfun(@(s)s(1),TE.ITIEnds(cellfun(@(s)~isempty(s),TE.ITIEnds))));
    df = df(~isnan(df));
    if all(df>0)   % monotonically increasing time stamps indicates that they are in abs. times
        for iT = 1:NumTrials
            TE.ITIEnds{iT} = TE.ITIEnds{iT} - oldTE.TrialStart(iT);   % make it rel. to TrialStart
        end
    end
    
    % Check signs
    if ~all(TE.TrialEnd(~isnan(TE.TrialEnd))>=0) || ...
            ~all(cell2mat(TE.ITIBegins')>=0) || ...
            ~all(cell2mat(TE.ITIEnds')>=0)
        error('fixTE:signError','Negative time stamps.')
    end
    
    % NoPokeITIBegins
    if isfield(TE,'NoPokeITIBegins')
        df = diff(cellfun(@(s)s(1),TE.NoPokeITIBegins(cellfun(@(s)~isempty(s),TE.NoPokeITIBegins))));
        df = df(~isnan(df));
        if all(df>0)   % monotonically increasing time stamps indicates that they are in abs. times
            for iT = 1:length(TE.NoPokeITIBegins)   % for some animal (18 to 23), last trial is missing
                TE.NoPokeITIBegins{iT} = TE.NoPokeITIBegins{iT} - oldTE.TrialStart(iT);   % make it rel. to TrialStart
            end
        end
        if ~all(cell2mat(TE.NoPokeITIBegins')>=0)   % check signs
            error('fixTE:signError','Negative time stamps.')
        end
    end
    
    % NoPokeITIEnds
    if isfield(TE,'NoPokeITIEnds')
        df = diff(cellfun(@(s)s(1),TE.NoPokeITIEnds(cellfun(@(s)~isempty(s),TE.NoPokeITIEnds))));
        df = df(~isnan(df));
        if all(df>0)   % monotonically increasing time stamps indicates that they are in abs. times
            for iT = 1:length(TE.NoPokeITIEnds)
                TE.NoPokeITIEnds{iT} = TE.NoPokeITIEnds{iT} - oldTE.TrialStart(iT);   % make it rel. to TrialStart
            end
        end
        if ~all(cell2mat(TE.NoPokeITIEnds')>=0)   % check signs
            error('fixTE:signError','Negative time stamps.')
        end
    end
    
    % LastITIBegins
    TE.LastITIBegins = nan(size(TE.ITIBegins));
    inx = cellfun(@(s)~isempty(s),TE.ITIBegins);   % in a session which is merged from 2, there is one empty in the middle
    TE.LastITIBegins(inx) = cellfun(@(s)s(end),TE.ITIBegins(inx));   % last ITI start
    
    % LastITIEnds
    TE.LastITIEnds = nan(size(TE.ITIEnds));
    inx = cellfun(@(s)~isempty(s),TE.ITIEnds);
    TE.LastITIEnds(inx) = cellfun(@(s)s(end),TE.ITIEnds(inx));   % last ITI ends
    
    % Save
%     save(fullfile(sessionpath,'TrialEvents.mat'),'-struct','TE')
end

% -------------------------------------------------------------------------
function TE2 = shortenTE(TE2,shinx)

% Eliminate behavioral trials
fnm = fieldnames(TE2);
for k = 1:length(fieldnames(TE2))
    TE2.(fnm{k}) = TE2.(fnm{k})(shinx);
end