function fr = firingrate(cellid,event_type,event,prewindow,postwindow,varargin)
%FIRINGRATE   Firing rate.
%   FR = FIRINGRATE(CELLID,EVENT_TYPE,EVENT,PREWINDOW,POSTWINDOW,VARARGIN)
%   calculates firing rate for a cell (CELLID) in windows aligned to a
%   reference event (EVENT).
%
%   Input arguments:
%       CELLID: cell ID in CellBase
%       EVENT_TYPE: 'stim' or 'trial' for stimulation and behavioral trial
%           events
%       EVENT: name of reference event
%       PREWINDOW: window boudaries before the event in seconds
%       POSTWINDOW: window boudaries after the event in seconds
%
%   Optional parameter value pairs (with default values):
%       'event_filter', 'none' - filter trials; see FILTERTRIALS for 
%           implemented filter types
%       'filterinput',[] - some filters require additional input; see
%           FILTERTRIALS for details
%
%   See also REGRESSION_ANALYSIS, FIRINGRATE_ANALYSIS and FIRINGRATE_CALL.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   24-Mar-2015

%   Edit log: BH, 3/24/15

% Default arguments
prs = inputParser;
addRequired(prs,'cellid',@iscellid)
addRequired(prs,'event_type',@ischar)   % event type ('stim' or 'trial')
addRequired(prs,'event',@(s)ischar(s)|...
    (iscellstr(s)&isequal(length(s),2))|...
    isa(s,'function_handle'))   % reference event
addRequired(prs,'prewindow',@(s)isnumeric(s)&isequal(length(s),2))  % time window relative to the event, in seconds
addRequired(prs,'postwindow',@(s)isnumeric(s)&isequal(length(s),2))  % time window relative to the event, in seconds
addParamValue(prs,'event_filter','none',@(s)ischar(s)|iscellstr(s))   % filter events based on properties
addParamValue(prs,'filterinput',[])   % some filters need additional input
parse(prs,cellid,event_type,event,prewindow,postwindow,varargin{:})
g = prs.Results;

% Load event structure
event_type = lower(event_type(1:4));
switch event_type
    case 'stim'
        
        % Load stimulus events
        try
            VE = loadcb(cellid,'StimEvents');   % load events
            VS = loadcb(cellid,'STIMSPIKES');   % load prealigned spikes
        catch ME
            disp('There was no stim protocol for this session.')
            error(ME.message)
        end
        
    case 'tria'
        
        % Load trial events
        try
            VE = loadcb(cellid,'TrialEvents');   % load events
            VS = loadcb(cellid,'EVENTSPIKES');   % load prealigned spikes
        catch ME
            disp('There was no behavioral protocol for this session.')
            error(ME.message)
        end
end

valid_trials = filterTrials(cellid,'event_type',event_type,'event',event,...
    'event_filter',g.event_filter,'filterinput',g.filterinput);   % filter trials
event_pos = findcellstr(VS.events(:,1),event);
spikes_stimulus = VS.event_stimes{event_pos}(valid_trials);   % spikes relative to ref. event

NumTrials = length(valid_trials);
fr = nan(NumTrials,2);   % firing rates
for iT = 1:NumTrials-1     % last trial may be incomplete
    lspikes = spikes_stimulus{iT};   % spikes relative to the event
    prespikes = lspikes(lspikes>prewindow(1)&lspikes<prewindow(2));   % spikes in the time window before the event
    postspikes = lspikes(lspikes>postwindow(1)&lspikes<postwindow(2));   % spikes in the time window after the event
    fr(iT,1) = length(prespikes) / diff(prewindow);   % firing rate in the time window before the event
    fr(iT,2) = length(postspikes) / diff(postwindow);   % firing rate in the time window after the event
end