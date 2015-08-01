function H = vippsth2_A1(cellid)
%VIPPSTH2_A1   Peri-stimulus time histogram.
%   H = VIPPSTH2_A1(CELLID) calculates PSTH aligned to light pulses
%   ('PulseOn') for the cell with the specified CellID. It returns the
%   handle of the resulting figure.
%
%   See also VIPISINFLUENCED3 and VIPPSTH_CALL.

% Input argument check
if nargin < 2
    win = [-0.25 0.65];  % time window for bin raster
    dt = 0.005;   % resolution of bin raster in s
    dsply = 1;   % supress display
end

% Set parameters and load CellBase variables
EventName1 = 'PulseOn';
ST = loadcb(cellid,'STIMSPIKES');   % load prealigned spikes for stimulation events
TE = loadcb(cellid,'StimEvents');
epoch_pos1 = findcellstr(ST.events(:,1),EventName1);
if epoch_pos1 == 0
    error('Epoch name not found');
end
stimes1 = ST.event_stimes{epoch_pos1};
time = win(1):dt:win(end);
valid_trials1 = find(~isnan(getfield(TE,EventName1)));

% Calculate bin rasters
spt1 = stimes2binraster(stimes1(valid_trials1),time,dt);

% Set input arguments for rater plot and PSTH
H = figure;
if dsply
    SEvent = 'PulseOn';
    FNum = 2;
    parts = 'all';
    sigma = 0.001;
    PSTHstd = 'on';
    
    % Plot raster plot and PSTH for 'BurstOn'
    set(gcf,'renderer','painters')   % temporaray change renderer because OpenGL locks the plot which result an error in legend layout handling
    viewcell2b(cellid,'TriggerName',EventName1,'SortEvent',SEvent,...
        'FigureNum',FNum,'eventtype','stim','window',win,'dt',dt,'sigma',sigma,'PSTHstd',PSTHstd,'Partitions',parts,...
        'EventMarkerWidth',0,'PlotZeroLine','on')
    pause(0.05)   % if reset the renderer two early, the same error occurs
    set(gcf,'renderer','opengl')   % reset renderer
end