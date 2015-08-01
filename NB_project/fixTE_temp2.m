function fixTE_temp2(cellids)
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
    [animal session] = cellid2tags(cellid);   % mouse ID and session ID
    sessionpath = fullfile(cellbase_path,animal,session);   % session full path   
    
    % LastITIBegins
    TE.LastITIBegins = nan(size(TE.ITIBegins));
    inx = cellfun(@(s)~isempty(s),TE.ITIBegins);   % in a session which is merged from 2, there is one empty in the middle
    TE.LastITIBegins(inx) = cellfun(@(s)s(end),TE.ITIBegins(inx));   % last ITI start
    
    % LastITIEnds
    TE.LastITIEnds = nan(size(TE.ITIEnds));
    inx = cellfun(@(s)~isempty(s),TE.ITIEnds);
    TE.LastITIEnds(inx) = cellfun(@(s)s(end),TE.ITIEnds(inx));   % last ITI ends
    
    % Save
    save(fullfile(sessionpath,'TrialEvents.mat'),'-struct','TE')
end

% -------------------------------------------------------------------------
function TE2 = shortenTE(TE2,shinx)

% Eliminate behavioral trials
fnm = fieldnames(TE2);
for k = 1:length(fieldnames(TE2))
    TE2.(fnm{k}) = TE2.(fnm{k})(shinx);
end