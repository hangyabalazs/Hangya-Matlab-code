%% QUICKANALYSIS2   Analysis of tetrode data.
%   QUICKANALYSIS2 is designed as an offline analysis tool for tetrode data
%   and behavior, which can be executed in an unsupervised manner on a
%   daily bases. It gives a quick overview of the experiment including
%   response profiles of clustered neurons, light-triggered PSTH and
%   psychometric plots of behavior. It relies on CellBase data handling
%   system.
%
%   See also INITCB.

% New to CellBase
add2cb = 0;     % 0, already in CellBase; 1, new to CellBase

% Stop if error
dbstop if error

% Animal, session
animalID = 'n028';
animalID2 = 'nb028';
sessionID = '120316b';
fullpth = [getpref('cellbase','datapath') '\' animalID '\' sessionID '\'];

% Directories
global DATAPATH
resdir = [DATAPATH 'OT\_response_profiles_temp\' animalID2 '\'];
if ~isdir(resdir)
    mkdir(resdir)
end
resdir2 = [DATAPATH 'OT\_behavior_temp\' animalID2 '\'];
if ~isdir(resdir2)
    mkdir(resdir2)
end

% Convert events file
nlxcsc2mat2(fullpth,'Channels','Events')

% Update CellBase
if add2cb
    addnewcells
end
cellids = findcell('rat',animalID,'session',sessionID);
disp(cellids)

% Light effects
if add2cb
    
    % Create stimulus events
    problem_stim_cellid = [];
    for iC = 1:length(cellids),
        cellid = cellids(iC);
        pathname = cellid2fnames(cellid,'Sess');
        try
            MakeStimEvents2(pathname,'BurstStartNttl',4)
        catch
            problem_stim_cellid = [problem_stim_cellid cellid];
        end
        try
            SE = load([pathname filesep 'StimEvents']);
            if isnan(SE.PulseOn),
                MakeStimEvents2(pathname,'BurstStartNttl',2)
            end
        catch
        end
    end
    
    % Prealign spikes to stimulus events
    problem_stim_cellid = [];
    for iC = 1:length(cellids)
        cellid = cellids(iC);
        try
            prealignSpikes(cellid,'FUNdefineEventsEpochs',@defineEventsEpochs_laserstim,'filetype','stim','ifsave',1,'ifappend',0)
        catch
            disp('Error in prealignSpikes.');
            problem_stim_cellid = [problem_stim_cellid cellid];
        end
    end
end

%% View light-triggered raster and PSTH
TrigEvent = 'BurstOn';
SEvent = 'BurstOn';
FNum = 2;
win = [-2 2.5];
% parts = 'all';
parts = '#BurstNPulse';
dt = 0.001;
sigma = 0.001;
PSTHstd = 'on';
ShEvent = {{'PulseOn','PulseOff','BurstOff'}};
ShEvColors = hsv(length(ShEvent{1}));
ShEvColors = mat2cell(ShEvColors,ones(size(ShEvColors,1),1),3);
for iCell = 1:length(cellids)
    cellid = cellids(iCell);
    H = figure;
    viewcell2b(cellid,'TriggerName',TrigEvent,'SortEvent',SEvent,'ShowEvents',ShEvent,'ShowEventsColors',{ShEvColors},...
        'FigureNum',FNum,'eventtype','stim','window',win,'dt',dt,'sigma',sigma,'PSTHstd',PSTHstd,'Partitions',parts,...
        'EventMarkerWidth',0,'PlotZeroLine','off')
    maximize_figure(H)
    
    cellidt = cellid{1};
    cellidt(cellidt=='.') = '_';
    fnm = [resdir cellidt '_LS.jpg'];   % save
    saveas(H,fnm)
    close(H)
end

%% View light-triggered raster and PSTH
TrigEvent = 'ZeroPulse';
SEvent = 'BurstOn';
FNum = 2;
win = [-2 2.5];
% parts = 'all';
parts = '#BurstNPulse';
dt = 0.001;
sigma = 0.001;
PSTHstd = 'on';
ShEvent = {{'PulseOn','PulseOff','BurstOff'}};
ShEvColors = hsv(length(ShEvent{1}));
ShEvColors = mat2cell(ShEvColors,ones(size(ShEvColors,1),1),3);
for iCell = 1:length(cellids)
    cellid = cellids(iCell);
    H = figure;
    viewcell2b(cellid,'TriggerName',TrigEvent,'SortEvent',SEvent,'ShowEvents',ShEvent,'ShowEventsColors',{ShEvColors},...
        'FigureNum',FNum,'eventtype','stim','window',win,'dt',dt,'sigma',sigma,'PSTHstd',PSTHstd,'Partitions',parts,...
        'EventMarkerWidth',0,'PlotZeroLine','off')
    maximize_figure(H)
    
    cellidt = cellid{1};
    cellidt(cellidt=='.') = '_';
    fnm = [resdir cellidt '_0LS.jpg'];   % save
    saveas(H,fnm)
    close(H)
end

%% Omission protocols
if sum(~isnan(SE.OmitPulse)) > 0
    
    % View light-triggered raster and PSTH
    TrigEvent = 'OmitPulse';
    SEvent = 'BurstOff';
    FNum = 2;
    win = [-1 1];
    % parts = 'all';
    parts = '#BurstNPulse';
    dt = 0.001;
    sigma = 0.001;
    PSTHstd = 'on';
    ShEvent = {{'PulseOn','PulseOff','BurstOff'}};
    ShEvColors = hsv(length(ShEvent{1}));
    ShEvColors = mat2cell(ShEvColors,ones(size(ShEvColors,1),1),3);
    for iCell = 1:length(cellids)
        cellid = cellids(iCell);
        H = figure;
        viewcell2b(cellid,'TriggerName',TrigEvent,'SortEvent',SEvent,'ShowEvents',ShEvent,'ShowEventsColors',{ShEvColors},...
            'FigureNum',FNum,'eventtype','stim','window',win,'dt',dt,'sigma',sigma,'PSTHstd',PSTHstd,'Partitions',parts,...
            'EventMarkerWidth',0,'PlotZeroLine','off')
        maximize_figure(H)
        
        cellidt = cellid{1};
        cellidt(cellidt=='.') = '_';
        fnm = [resdir cellidt '_OSR.jpg'];   % save
        saveas(H,fnm)
        close(H)
    end
end

%% Cluster quality
if add2cb
    BatchSessionClust(fullpth)
end

%% Separate protocols

% Trigger event
trig_event = 'ZeroPulse';

% Load stimulus events
SE = loadcb(cellid,'StimEvents');

% Addition to allow partition based on ProtocolID
for k=1:length(SE.ProtocolID)
    SE.ProtocolID2(k) = str2num(SE.ProtocolID{k}(2:end));
end

% Addition to allow ShowEvents for PulseOn
boinx = find(~isnan(SE.(trig_event)));
pulseon_cell = cell(1,length(SE.(trig_event)));
if isequal(trig_event,'ZeroPulse')
    pulsediff = diff(SE.PulseOn);
else
    pulsediff = 0;
end
for k = 1:length(boinx)-1
    pinx = boinx(k):(boinx(k+1)-1);
    pulseon_cell{boinx(k)} = SE.PulseOn(pinx)' - pulsediff(pinx(1));
end
pulseon_cell{end} = SE.PulseOn(boinx(end):end);
SE.PulseOnCell = pulseon_cell;

%%

% eventtype stim
% TrigEvent='BurstOn';
TrigEvent = 'ZeroPulse';
% TrigEvent='OmitPulse';
SEvent='BurstOff';
FNum=2;
win=[-0.2 2.5];
parts='all';
% parts='#BurstNPulse';
parts='#ProtocolID2';
dt=0.001;
sigma=0.001;
PSTHstd='on';
ShEvent={{'PulseOnCell'}};
ShEvColors=hsv(length(ShEvent{1}));
ShEvColors=mat2cell(ShEvColors,ones(size(ShEvColors,1),1),3);
for iCell=1:length(cellids),
    cellid=cellids(iCell);
%     cellid=allcells(iCell);
%     try
        figure
        viewcell2b(cellid,'TriggerName',TrigEvent,'SortEvent',SEvent,'ShowEvents',ShEvent,'ShowEventsColors',{ShEvColors},...
        'FigureNum',FNum,'eventtype','stim','window',win,'dt',dt,'sigma',sigma,'PSTHstd',PSTHstd,'Partitions',parts,...
        'EventMarkerWidth',0,'PlotZeroLine','off')
%         pause(1)
%     catch
%     end
end