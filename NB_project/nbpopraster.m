function nbpopraster(cellids,varargin)
%NBPOPRASTER   Population raster.
%   NBPOPRASTER(CELLIDS) calculates population raster plot for a set of
%   cells (CELLIDS - see CellBase documentation for details on cell ID
%   form). Rasters of single cells of the same trial are aligned to the
%   same event.
%
%   Default behavior of NBPOPRASTER can be modified by using a set of
%   paramter-value pairs as optional input parameters. The following
%   parameters are implemented (with default values):
%   	'event_type', 'stim' - type of reference event, 'stim' or 'trial'
%           for stimulus and trial events, respectively
%       'TriggerEvent', 'BurstOn' - reference event; default is 'BurstOn'
%       'ShowEvents', 'PulseOnCell' - events to mark on the resulting
%           raster; default behavior marks all stimulus presentations
%       'parts', 'all' - partitioning the set of trials; input to
%           PARTITION_TRIALS, see details therein (default, no
%           partitioning)
%       'eventindex', 1 - index of concrete, single reference event within
%           the event stream
%   	'window', [-0.5 3] - time window for the raster, time relative to 0 
%           (reference event) in seconds
%   	'dt', 0.001 - time resolution in seconds
%
%   See also ULTIMATE_PSTH and NBPOPPSTH.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   24-Oct-2012

%   Edit log: BH 10/24/12

% Default arguments
prs = inputParser;
addRequired(prs,'cellids',@(s)ischar(s)|iscellstr(s))   % cellID
addParamValue(prs,'event_type','stim',@ischar)   % event type ('stim' or 'trial')
addParamValue(prs,'TriggerEvent','BurstOn',@(s)ischar(s))   % default reference event: 'BurstOn'
addParamValue(prs,'ShowEvents',{'PulseOnCell'},@(s)iscellstr(s))   % default: show light pulses
addParamValue(prs,'eventindex',1,@(s)isnumeric(s))   % selection for a particular event
addParamValue(prs,'window',[-0.5 3],@(s)isnumeric(s)&isequal(length(s),2))  % time window relative to the event, in seconds
addParamValue(prs,'dt',0.001,@isnumeric)   % time resolution, in seconds
parse(prs,cellids,varargin{:})
g = prs.Results;
if ischar(cellids)   % convert 'cellids' to cell
    cellids = {cellids};
end

% Load events
try
    TE = loadcb(cellids{1},'StimEvents');   % load events
catch ME
    disp('There was no behavioral protocol for ths session.')
    error(ME.message)
end

% Load spikes
NumCells = length(cellids);
for iC = 1:NumCells
    ST = loadcb(cellids{iC},'STIMSPIKES');   % load prealigned spikes
    
    stim_pos = findcellstr(ST.events(:,1),g.TriggerEvent);   % tone onset
    spikes_stimulus = ST.event_stimes{stim_pos};   % spikes relative to burst onset
    stimes(iC) = spikes_stimulus(g.eventindex);   % select one burst 
end

% Addition to allow ShowEvents for PulseOn
time = g.window(1):g.dt:g.window(2);   % time base array
boinx = find(~isnan(TE.(g.TriggerEvent)));
pulseon_cell = cell(1,length(TE.(g.TriggerEvent)));
pulsediff = zeros(length(TE.PulseOn));
for k = 1:length(boinx)-1
    pinx = boinx(k):(boinx(k+1)-1);
    pulseon_cell{boinx(k)} = TE.PulseOn(pinx)' - pulsediff(pinx(1));
end
pulseon_cell{end} = TE.PulseOn(boinx(end):end);
TE.PulseOnCell = pulseon_cell;

% Show light pulses
EventTimes = trialevents2relativetime(TE,g.TriggerEvent,g.ShowEvents);
et = EventTimes{g.eventindex}';
det = 0.001;
et(et<time(1)) = NaN;
et(et>time(end)) = NaN;
H = figure;
hold on
set(gca,'Color',[0.15 0.15 0.15])
p = patch([et; et; et+det; et+det],repmat([NumCells+0.5; 0.5; 0.5; NumCells+0.5],1,size(et,2)),'k');
set(p,'FaceColor',[0 153 255]/255,'EdgeColor',[0 153 255]/255,'FaceAlpha',0.5);

% Raster plot
for iC = 1:NumCells
    spikes = stimes{iC};
    ok_spikes = spikes(spikes>g.window(1)&spikes<=g.window(2));
    ok_spikes = ok_spikes(:)';
    if ~isempty(ok_spikes)
        line([ok_spikes; ok_spikes],[(iC-0.28)*ones(1,length(ok_spikes)); (iC+0.28)*ones(1,length(ok_spikes))],...
            'Color','w','Linewidth',5);
    end
end

% Store parameters with figure
setappdata(H,'cellids',cellids)
setappdata(H,'parameters',g)
keyboard