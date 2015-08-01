function [R L J B M lim1 lim2] = reliability_latency_jitter(cellid,varargin)
%RELIABILITY_LATENCY_JITTER   Reliability, latency and jitter of spiking.
%   [R L J] = RELIABILITY_LATENCY_JITTER(CELLID) calculates the latency of
%   spiking of a cell (CELLID) after an event (default, 'PulseOn') (L). It
%   also returns the standard deviation of spike times after finding the
%   period of stimulated spikes (see FINDSTIMPERIOD for details) (J).
%   Reliability of evoking spikes (proportion of events followed by evoked
%   spikes) is returned in R.
%
%   [R L J] = RELIABILITY_LATENCY_JITTER(CELLID,'EVENT',EV) uses EV (e.g.
%   'BurstOn') as reference event instead of 'PulseOn'.
%
%   [R L J] = RELIABILITY_LATENCY_JITTER(CELLID,'EVENT',EV,...
%       'EVENT_FILTER',FILTER,'FILTERINPUT',EI)
%   filters the trials using FILTEREVENT function (see details on possible
%   filters there). Some filters need additional input that can be provided
%   using the 'FILTERINPUT' parameter.
%
%   [R L J B M] = RELIABILITY_LATENCY_JITTER(CELLID,...) also returns
%   baseline firing rate (B) and maximal firing rate in the test window
%   (M).
%
%   [R L J B M A1 A2] = RELIABILITY_LATENCY_JITTER(CELLID,...) returns the
%   start and end point of detected stimulation period (see
%   FINDSTIMPERIOD).
%
%   Examples:
%   [E L J] = reliability_latency_jitter(cellid,'event','BurstOn');
%
%   fi = struct('BurstNPulse',20);
%   [E I J] = reliability_latency_jitter(cellid,...
%       'event_filter','BurstNPulse_maxPower','filterinput',fi);
%
%   [E_burston L_burston J_burston B_burston M_burston Astart Aend] = ...
%        reliability_latency_jitter(cellid,'event','PulseOn');
%
%   See also FINDSTIMPERIOD, EXTRACTSEGSPIKES and EFFICIENCY.

% Input arguments
prs = inputParser;
addRequired(prs,'cellid',@iscellid)
addParamValue(prs,'event','PulseOn',@ischar)   % reference event
addParamValue(prs,'event_filter','none',@ischar)   % event filter
addParamValue(prs,'filterinput',[])   % some filters need additional input
parse(prs,cellid,varargin{:})
g = prs.Results;

% Filter events
valid_trials = filterTrials(cellid,'event_type','stim','event',g.event,...
    'event_filter',g.event_filter,'filterinput',g.filterinput);

% Find stimulation peak (latency)
[lim1 lim2 L at B M] = findStimPeriod(cellid,'event',g.event,...
    'valid_trials',valid_trials);   %#ok<ASGLU> % find stimulated period and peak time

% Standard deviation of evoked spikes (jitter)
tsegs_evoked = rel2abstimes(cellid,[lim1 lim2],'stim',g.event,'valid_trials',valid_trials);   % convert period to epochs relative to event
selts_evoked = extractSegSpikes(cellid,tsegs_evoked);   % find evoked spikes
relevokedtimes = abs2reltimes(cellid,selts_evoked,'stim',g.event);  % convert absolute spike times to times rel. to the event
J = std(relevokedtimes);

% Reliability
R = efficiency(cellid,selts_evoked,'stim',g.event,'valid_trials',valid_trials);   % efficiency of evoking spikes