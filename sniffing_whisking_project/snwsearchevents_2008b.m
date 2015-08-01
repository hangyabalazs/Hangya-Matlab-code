function [matchedEvents matchedTS matchedString] = snwsearchevents_2008b(str)
%SNWSEARCHEVENTS   Perform search on Neuralynx events.
%   SNWSEARCHEVENTS supports two search modes.
%   (1) Search for a character array in Neuralynx event strings.
%   ME = SNWSEARCHEVENTS(STR) looks for the string STR (or 'eating' if no 
%   input argument is specified) within Neuralynx event strings. Path names
%   of sessions with at least one match are returned in ME. SNWSEARCHEVENTS
%   also supports cell array of strings as input argument. A match is
%   detected in case any of the strings in the cell array is contained in
%   any of the event strings. Note that the search is not case sensitive.
%
%   Example:
%   matchedEvents = snwsearchevents({'eating' 'food'});
%
%   (2) Search for a specific event ID in Neuralynx event IDs. 
%   ME = SNWSEARCHEVENTS(K) looks for IDs with a value of K among Neuralynx
%   event IDs. Path names of sessions with at least one match are returned
%   in ME. SNWSEARCHEVENTS also supports an array of ID values as input
%   argument. A match is detected in case any ID values in the array is 
%   matched with any of the event IDs.
%
%   Example:
%   matchedEvents = snwsearchevents([1 4]);
%
%   [ME MTS] = SNWSEARCHEVENTS(K) also returns corresponding timestamps of
%   matches (MTS) in a cell array. First index corresponds to the session
%   name in ME, whereas the second index corresponds to the entity in the
%   input argument list.
%
%   Example:
%   [matchedEvents matchedTS] = snwsearchevents([1 4]);
%
%   To reach the timestamps for the third matching session, i.e.
%   matchedEvents{3} with event ID value of 4, use cell array indexing as
%   follows; matchedTS{3,2}.
%
%   [ME MTS MS] = SNWSEARCHEVENTS(K) returns corresponding event strings of
%   matches (MS) in a cell array. Indexing is similar to that of MTS (see
%   above).
%
%   Example:
%   [matchedEvents matchedTS matchedString] = snwsearchevents([1 4]);
%
%   See also NLXEVENTS.

% Input argument check
error(nargchk(0,1,nargin))
if nargin < 1
    str = {'eating'};     % set default string
    type = 'EventString';
else
    if isnumeric(str)     % two types of search (for IDs or for strings)
        type = 'EventID';
    else
        type = 'EventString';
        if ~iscell(str)     % convert input to cell
            tmp = str;
            str = {tmp};
        end
    end
end
ls = length(str);

% Add Nlx converter to path
addpath('C:\MATLAB_R2010a\work\Nlx_converter')

% Import events
global DATADIR
inpdir = [DATADIR 'HSW' filesep];
dr = dir(inpdir);
dr = getsubdir(dr);     % select subdirectories
matchedEvents = {};     % initialize output
matchedTS = cell(0,ls);
matchedString = cell(0,ls);
for k = 1:length(dr)    % animal loop
    locdir = [inpdir dr(k).name filesep];
    dr2 = dir(locdir);
    dr2 = getsubdir(dr2);     % select subdirectories
    for k2 = 1:length(dr2)    % session loop
        fn = [locdir dr2(k2).name filesep];
        events = loadevents(fn);
        
        % Search NlxEvents object array
        switch type
            case 'EventString'
                evstrs = {events.EventString};
                evstrs = cellfun(@(s)lower(s),evstrs,'UniformOutput',0);  % case-insensitive search
                str = cellfun(@(s)lower(s),str,'UniformOutput',0);
                sfd = cellfun(@(s)strfind(evstrs,s),str,'UniformOutput',0);
                sfd2 = cellfun(@(s)cell2mat(s),sfd,'UniformOutput',0);
                if ~isempty(cell2mat(sfd2))
                    matchedEvents{end+1} = fn;
                    sfd2i = arrayfun(@(t)cell2mat(cellfun(@(s)isempty(s),sfd{t},'UniformOutput',0)),(1:ls),'UniformOutput',0);
                    pts = cellfun(@(s)find(~s),sfd2i,'UniformOutput',0);
                    matchedTS(end+1,:) = arrayfun(@(t)eventsref(events,pts{t},'TimeStamp'),(1:ls),'UniformOutput',0);
                    matchedString(end+1,:) = arrayfun(@(t)eventsrefs(events,pts{t},'EventString'),(1:ls),'UniformOutput',0);
                end
            case 'EventID'
                evids = [events.EventID];
                sfd = arrayfun(@(s)isempty(find(evids==s,1)),str);
                if any(~sfd)
                    matchedEvents{end+1} = fn;
                    pts = arrayfun(@(s)find(evids==s),str,'UniformOutput',0);
                    matchedTS(end+1,:) = arrayfun(@(t)eventsref(events,pts{t},'TimeStamp'),(1:ls),'UniformOutput',0);
                    matchedString(end+1,:) = arrayfun(@(t)eventsrefs(events,pts{t},'EventString'),(1:ls),'UniformOutput',0);
                end
        end
    end    % session loop
end    % animal loop
matchedEvents = matchedEvents';     % return as column

% -------------------------------------------------------------------------
function dr = getsubdir(dr)
idr = cell2mat({dr.isdir});
dr = dr(idr);       % select directories
dr = dr(3:end);     % select subdirectories
    
% -------------------------------------------------------------------------
function events = loadevents(fn)

% Load Nlx events and create objects
[EventTimeStamps, EventIDs, Nttls, Extras, EventStrings NlxHeader] = ...
    Nlx2MatEV([fn  'Events.nev'],[1 1 1 1 1],1,1,1);         % load Neuralynx events
numev = length(EventIDs);   % number of events
events = NlxEvents.empty(0,numev);
for k = 1:numev
    events(k) = NlxEvents(EventIDs(k),EventTimeStamps(k)/1000000,Nttls(k),EventStrings{k});
end

% -------------------------------------------------------------------------
function tt = eventsref(events,inx,fldname)

% Reference for numeric values in events object for maintaining compatibility
if isempty(inx)
    tt = [];
else
    tt = [events(inx).(fldname)];
end

% -------------------------------------------------------------------------
function tt = eventsrefs(events,inx,fldname)

% Reference for strings in events object for maintaining compatibility
if isempty(inx)
    tt = [];
else
    tt = {events(inx).(fldname)};
end