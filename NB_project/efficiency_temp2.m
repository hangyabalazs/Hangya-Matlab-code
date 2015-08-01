function [E spnodist stimes1st] = efficiency(cellid,stimes,event_type,event,varargin)
%EFFICIENCY   Efficiency of evoking spikes.
%   E = EFFICIENCY(CELLID,STIMES,EVENT_TYPE,EVENT) calculates the
%   proportion of events (EVENT) that were followed by a spike (STIMES) of
%   a given cell (CELLID). EVENT_TYPE should determine whether EVENT is a
%   field of stimulus events ('stim') or trial events ('trial').
%
%   E = EFFICIENCY(CELLID,STIMES,EVENT_TYPE,EVENT,'VALID_TRIALS',VT)
%   filters the trials according to the 'VALID_TRIALS', TR parameter value
%   pair.
%
%   [E SPNO] = EFFICIENCY(CELLID,STIMES,EVENT_TYPE,EVENT,'VALID_TRIALS',VT)
%   also returns the number of spikes detected for all valid trials (SPNO).
%
%   [E SPNO T1] = EFFICIENCY(CELLID,STIMES,EVENT_TYPE,EVENT,'VALID_TRIALS',VT)
%   returns a vector of spike times including only the first spike in each
%   session (T1).
%
%   See also RELIABILITY_LATENCY_JITTER.

% Input arguments
prs = inputParser;
addRequired(prs,'cellid',@iscellid)
addRequired(prs,'stimes',@isnumeric)  % time stamps
addRequired(prs,'event_type',@ischar)   % event type ('stim' or 'trial')
addRequired(prs,'event',@ischar)   % reference event
addParamValue(prs,'valid_trials','all',@(s)isnumeric(s)||islogical(s))  % valid trials - use all trials by default
parse(prs,cellid,stimes,event_type,event,varargin{:})
g = prs.Results;

% Load event structure
event_type = lower(event_type(1:4));
switch event_type
    case 'stim'
        
        % Load stimulus events
        try
            VE = loadcb(cellid,'StimEvents');
        catch ME
            disp('There was no stim protocol for ths session.')
            error(ME.message)
        end
        
    case 'tria'
        
        % Load trial events
        try
            VE = loadcb(cellid,'TrialEvents');
        catch ME
            disp('There was no behavioral protocol for ths session.')
            error(ME.message)
        end
        
    otherwise
        error('Input argument ''event_type'' should be either ''stim'' or ''trial''.')
end

% Valid trials
valid_trials = parseValidTrials(VE,event,g.valid_trials);

% Find the corresponding event for each spike
tno = length(stimes);
[evinx stimes1st] = deal(nan(1,tno));
stimes1st(1) = stimes(1);
next = 2;
for k = 1:tno
    if isequal(event_type,'tria') && ~isequal(event,'TrialStart')
        pr = stimes(k) - (VE.(event)(valid_trials) + VE.TrialStart(valid_trials));     % trial events other than TrialStart are stored relative to TrialStart
    else
        pr = stimes(k) - VE.(event)(valid_trials);
    end
    evinx(k) = find(pr>0,1,'last');
    if k > 1 && ~isequal(evinx(k),evinx(k-1))
        stimes1st(next) = stimes(k);   % spike times including first spikes only in each trial
        next = next + 1;
    end
end
stimes1st(next:end) = [];
spnodist = arrayfun(@(s)sum(evinx==s),1:length(valid_trials));   % number of spikes for all valid trials
evinx = unique(evinx);    % all events which were followed by a spike

% Efficiency
allev = VE.(event)(valid_trials);
allev = allev(~isnan(allev));   % all events
E = length(evinx) / length(allev);  % efficiency of evoking spikes