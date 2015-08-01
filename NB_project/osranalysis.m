function osranalysis
%QUICKANALYSIS2   Analysis of tetrode data.
%   QUICKANALYSIS2 is designed as an offline analysis tool for tetrode data
%   and behavior, which can be executed in an unsupervised manner on a
%   daily bases. It gives a quick overview of the experiment including
%   response profiles of clustered neurons, light-triggered PSTH and
%   psychometric plots of behavior. It relies on CellBase data handling
%   system.
%
%   See also INITCB.

% Animal, session
animalID = 'n023';
animalID2 = 'nb023';
sessionID = '120104b';
fullpth = [getpref('cellbase','datapath') '\' animalID '\' sessionID '\'];
% cellids = findcell('rat',animalID,'session',sessionID);
cellids = {'n023_120104b_4.1'};
disp(cellids)

% Directories
% global DATAPATH
% resdir = [DATAPATH 'OT\_response_profiles\' animalID2 '\'];
% if ~isdir(resdir)
%     mkdir(resdir)
% end
% resdir2 = [DATAPATH 'OT\_behavior\' animalID2 '\'];
% if ~isdir(resdir2)
%     mkdir(resdir2)
% end

% Light effects
% View light-triggered raster and PSTH
TrigEvent = 'BurstOn';
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
    
%     cellidt = cellid{1};
%     cellidt(cellidt=='.') = '_';
%     fnm = [resdir cellidt '_LS.jpg'];   % save
%     saveas(H,fnm)
%     close(H)
    
    t2l = [];
    t2o = [];
    prd = [];
    
    
    keyboard
    
    H2 = figure;
    axes
    
    
    x = get(gco,'XData');
    y = get(gco,'YData');
%     figure;
%     plot(x,smooth(y,'linear',5))
    [vdisc thr] = b_udisc(smooth(y,'linear',5));
    osr_pulse = linterp(1:length(x),x,vdisc(end));
    num_pulses = input('Number of pulses?')
    freq = num_pulses / 2;
    stimulus_period = 1 / freq;
    omitted_pulse = num_pulses * stimulus_period;
    last_pulse = (num_pulses - 1) * stimulus_period;
    time2last = osr_pulse - last_pulse;
    time2omitted = osr_pulse - omitted_pulse;
    
    prd = [prd stimulus_period];
    t2l = [t2l time2last];
    t2o = [t2o time2omitted];
    
    figure;plot(prd,t2l-t2o)
    figure;plot(prd,t2l,'bo-')
    hold on;plot(prd,t2o,'ko-')
end