function [psth spsth spsth_se] = ultimate_psth(cellid,event_type,event,window,varargin)
%ULTIMATE_PSTH   Peri-stimulus time histogram.
%   [PSTH SPSTH SPSTH_SE] = ULTIMATE_PSTH(CELLID,EVENT_TYPE,EVENT,WINDOW,VARARGIN)
%   calculates peri-stimulus time histogram (PSTH) for the cell passed in
%   CELLID. Smoothed PSTH (SPSTH) and SE of smoothing (SPSTH_SE) are also
%   returned.
%
%   Mandatory input arguments:
%       CELLID: defines the cell (see CellBase documentation)
%       EVENT: the event to which the PSTH is aligned
%       EVENT_TYPE: the type of event, 'stim' or 'trial'
%       WINDOW: window for calculation relative to the event in seconds
%
%   Default behavior of ULTIMATE_PSTH can be modified by using a set of
%   paramter-value pairs as optional input parameters. The following
%   parameters are implemented (with default values):
%   	'dt', 0.001 - time resolution in seconds
%       'sigma', 0.02 - smoothing kernel for the smoothed PSTH, in seconds
%       'margin',[-0.01 0.01] margins for PSTH calculation to get rid of 
%           edge effect due to smoothing
%       'event_filter', 'none' - filter light-stimulation trials; see
%           FILTERTRIALS for implemented filter types
%       'filterinput',[] - some filters require additional input; see
%           FILTERTRIALS for details
%       'maxtrialno', 5000 - maximal number of trials included; if ther are
%           more valid trials, they are randomly down-sampled
%       'parts', 'all' - partitioning the set of trials; input to
%           PARTITION_TRIALS, see details therein (default, no
%           partitioning)
%       'isadaptive, true - if false, classic PSTH algorithm is applied; if
%           true, adaptive PSTH is calculated (see APSTH2)
%
%   See also BINRASTER2PSTH, BINRASTER2APSTH, APSTH2, VIEWCELL2B,
%   PARTITION_TRIALS and FILTERTRIALS.

%   Edit log: BH 7/5/12

% Default arguments
prs = inputParser;
addRequired(prs,'cellid',@iscellid)
addRequired(prs,'event_type',@ischar)   % event type ('stim' or 'trial')
addRequired(prs,'event',@ischar)   % default reference event: 'LeftPortIn'
addRequired(prs,'window',@(s)isnumeric(s)&isequal(length(s),2))  % time window relative to the event, in seconds
addParamValue(prs,'event_filter','none',@ischar)   % filter events based on properties
addParamValue(prs,'filterinput',[])   % some filters need additional input
addParamValue(prs,'maxtrialno',5000)   % downsample events if more than 'maxtrialno'
addParamValue(prs,'dt',0.001,@isnumeric)   % time resolution of the binraster, in seconds
addParamValue(prs,'sigma',0.02,@isnumeric)     % smoothing kernel for the smoothed PSTH
addParamValue(prs,'margin',[-0.1 0.1])  % margins for PSTH calculation to get rid of edge effect due to smoothing
addParamValue(prs,'parts','all')   % partition trials
addParamValue(prs,'isadaptive',true,@(s)islogical(s)|ismember(s,[0 1]))   % use adaptive PSTH algorithm
addParamValue(prs,'display',false,@(s)islogical(s)|ismember(s,[0 1]))   % control displaying rasters and PSTHs
parse(prs,cellid,event_type,event,window,varargin{:})
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
            disp('There was no stim protocol for ths session.')
            error(ME.message)
        end
        
    case 'tria'
        
        % Load trial events
        try
            VE = loadcb(cellid,'TrialEvents');   % load events
            VS = loadcb(cellid,'EVENTSPIKES');   % load prealigned spikes
        catch ME
            disp('There was no behavioral protocol for ths session.')
            error(ME.message)
        end
        
    otherwise
        error('Input argument ''event_type'' should be either ''stim'' or ''trial''.')
end

% Set parameters and load CellBase variables
event_pos = findcellstr(VS.events(:,1),event);
if event_pos == 0
    error('Event name not found');
end
stimes = VS.event_stimes{event_pos};
time = (window(1)+g.margin(1)):g.dt:(window(2)+g.margin(2));
valid_trials = filterTrials(cellid,'event_type',event_type,'event',event,...
    'event_filter',g.event_filter,'filterinput',g.filterinput);

% Calculate bin rasters
spt = stimes2binraster(stimes,time,g.dt);

% Partition trials
[COMPTRIALS, TAGS] = partition_trials(VE,g.parts);

% PSTH
if g.isadaptive
    [psth, spsth, spsth_se] = binraster2apsth_temp(spt,g.dt,g.sigma,COMPTRIALS,valid_trials);
else
    [psth, spsth, spsth_se] = binraster2psth(spt,g.dt,g.sigma,COMPTRIALS,valid_trials);
end
stm = abs(window(1)+g.margin(1)) / g.dt;   % zero point
inx = (stm+1+window(1)/g.dt):(stm+1+window(2)/g.dt);     % indices for cutting margin
psth = psth(inx);
spsth = spsth(inx);
spsth_se = spsth_se(inx);