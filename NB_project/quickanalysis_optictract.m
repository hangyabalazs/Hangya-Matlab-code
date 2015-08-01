function quickanalysis_optictract
%QUICKANALYSIS_OPTICTRACT   Analysis of tetrode data.
%   QUICKANALYSIS_OPTICTRACT is designed as an offline analysis tool for
%   tetrode data and behavior, which can be executed in an unsupervised
%   manner on a daily bases. It gives a quick overview of the experiment
%   including response profiles of clustered neurons, light-triggered PSTH
%   and psychometric plots of behavior. It relies on CellBase data handling
%   system.
%
%   See also INITCB.

% New to CellBase
add2cb = 0;     % 0, already in CellBase; 1, new to CellBase

% Stop if error
dbstop if error

% Animal, session
animalID = 'n023';
animalID2 = 'nb023';
sessionID = '111229d';
fullpth = [getpref('cellbase','datapath') '\' animalID '\' sessionID '\'];

% Directories
global DATAPATH
resdir = [DATAPATH 'OT\_response_profiles2\' animalID2 '\'];
if ~isdir(resdir)
    mkdir(resdir)
end
resdir2 = [DATAPATH 'OT\_behavior\' animalID2 '\'];
if ~isdir(resdir2)
    mkdir(resdir2)
end

% Convert events file
if add2cb
    nlxcsc2mat2(fullpth,'Channels','Events')
end

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

% View light-triggered raster and PSTH
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

% View light-triggered raster and PSTH
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

% Omission protocols
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

% Cluster quality
if add2cb
    BatchSessionClust(fullpth)
end