function quickanalysis_optictract2(animalNO,sessionID,add2cb)
%QUICKANALYSIS_OPTICTRACT2   Analysis of tetrode data.
%   QUICKANALYSIS_OPTICTRACT2 is designed as an offline analysis tool for
%   tetrode data and behavior, which can be executed in an unsupervised
%   manner on a daily bases. It gives a quick overview of the experiment
%   including response profiles of clustered neurons, light-triggered PSTH
%   and psychometric plots of behavior. It relies on CellBase data handling
%   system.
%
%   See also INITCB.

% Stop if error
dbstop if error

% Input argument check
error(nargchk(0,3,nargin))

% New to CellBase
if nargin < 3
    add2cb = 0;     % 0, already in CellBase; 1, new to CellBase
end

% Animal, session
if nargin < 2
    sessionID = '120301b';
end
if nargin < 1
    animalID2 = 'nb023';
    animalID = 'n023';
else
    animalID2 = ['nb0' num2str(animalNO)];
    animalID = ['n0' num2str(animalNO)];
end
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

% Separate protocols
% Load stimulus events
SE = loadcb(cellids(1),'StimEvents');

% Trigger event
if isfield(SE,'ZeroPulse')
    trig_event = 'ZeroPulse';
else
    trig_event = 'BurstOn';
end

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
    pulsediff = zeros(length(SE.PulseOn));
end
for k = 1:length(boinx)-1
    pinx = boinx(k):(boinx(k+1)-1);
    pulseon_cell{boinx(k)} = SE.PulseOn(pinx)' - pulsediff(pinx(1));
end
pulseon_cell{end} = SE.PulseOn(boinx(end):end);
SE.PulseOnCell = pulseon_cell;

% Save stimulus events
fname = cellid2fnames(cellids(1),'StimEvents');
save(fname,'-struct','SE')

% Plot
TrigEvent = trig_event;
% TrigEvent='BurstOn';
% TrigEvent = 'ZeroPulse';
% TrigEvent='OmitPulse';
SEvent = 'BurstOff';
FNum = 2;
win=[-0.2 2.5];
parts='all';
parts = '#BurstNPulse';
% parts='#ProtocolID2';
dt = 0.001;
sigma = 0.001;
PSTHstd = 'on';
ShEvent = {{'PulseOnCell'}};
ShEvColors = hsv(length(ShEvent{1}));
ShEvColors = mat2cell(ShEvColors,ones(size(ShEvColors,1),1),3);
for iCell = 1:length(cellids),
    cellid = cellids(iCell);
    H = figure;
    viewcell2b(cellid,'TriggerName',TrigEvent,'SortEvent',SEvent,'ShowEvents',ShEvent,'ShowEventsColors',{ShEvColors},...
        'FigureNum',FNum,'eventtype','stim','window',win,'dt',dt,'sigma',sigma,'PSTHstd',PSTHstd,'Partitions',parts,...
        'EventMarkerWidth',0.008,'PlotZeroLine','off')
    maximize_figure(H)
        
    cellidt = cellid{1};
    cellidt(cellidt=='.') = '_';
    fnm = [resdir cellidt '_LS.jpg'];   % save
    saveas(H,fnm)
    close(H)
end

% Cluster quality
if add2cb
    BatchSessionClust(fullpth)
end