function quickanalysis2ro
%QUICKANALYSIS2RO   Analysis of tetrode data.
%   QUICKANALYSIS2RO is designed as an offline analysis tool for tetrode
%   data and behavior, which can be executed in an unsupervised manner on a
%   daily bases. It gives a quick overview of the experiment including
%   response profiles of clustered neurons, light-triggered PSTH and
%   psychometric plots of behavior. It relies on CellBase data handling
%   system.
%
%   QUICKANALYSIS2RO calculates response profiles for session already
%   incorporated in CellBase.
%
%   See also INITCB.

% Stop if error
dbstop if error

% Animal, session
animalID = 'n013';
animalID2 = 'nb013';
sessionID = '110630a';
fullpth = [getpref('cellbase','datapath') '\' animalID '\' sessionID '\'];

% Directories
global DATAPATH
resdir = [DATAPATH 'NB\_response_profiles\' animalID2 '\'];
if ~isdir(resdir)
    mkdir(resdir)
end
resdir2 = [DATAPATH 'NB\_behavior\' animalID2 '\'];
if ~isdir(resdir2)
    mkdir(resdir2)
end

% Get cell IDs
cellids = findcell('rat',animalID,'session',sessionID);

% Response profiles
% Is predictive?
for k = 1:length(cellids)
    H = figure;
%     viewcell2b(cellids(k),'TriggerName','StimulusOn','SortEvent','Stimulu
%     sOff','eventtype','behav','ShowEvents',{{'StimulusOff'}},'Partitions','#ResponseType','window',[-5 5])
    viewcell2b(cellids(k),'TriggerName','StimulusOn','SortEvent','TrialStart','eventtype','behav','ShowEvents',{{'StimulusOff'}},'Partitions','#ResponseType','window',[-5 5])
    maximize_figure(H)
    
    cellidt = cellids{k};
    cellidt(cellidt=='.') = '_';
    fnm = [resdir cellidt '_IPD.jpg'];   % save
    saveas(H,fnm)
    close(H)
end

% Hit & FA
for k = 1:length(cellids)
    H = figure;
    viewcell2b(cellids(k),'TriggerName','LeftPortIn','SortEvent','StimulusOn','eventtype','behav','ShowEvents',{{'StimulusOn'}},'Partitions','#ResponseType','window',[-5 5])
    maximize_figure(H)
    
    cellidt = cellids{k};
    cellidt(cellidt=='.') = '_';
    fnm = [resdir cellidt '_HF.jpg'];   % save
    saveas(H,fnm)
    close(H)
end

% Does it depend on stim intensity?
for k = 1:length(cellids)
    H = figure;
    viewcell2b(cellids(k),'TriggerName','StimulusOn','SortEvent','StimulusOff','eventtype','behav','ShowEvents',{{'StimulusOff'}},'Partitions','#StimulusDuration','window',[-5 5])
    maximize_figure(H)
    
    cellidt = cellids{k};
    cellidt(cellidt=='.') = '_';
    fnm = [resdir cellidt '_SI.jpg'];   % save
    saveas(H,fnm)
    close(H)
end

% Lickraster
% for k = 1:length(cellids)
%     H = figure;
%     viewcell2b(cellids(k),'TriggerName','LickIn','SortEvent','StimulusOff','eventtype','behav','ShowEvents',{{'StimulusOff'}},'Partitions','#ResponseType','window',[-5 5])
%     maximize_figure(H)
%     
%     fnm = [resdir cellid '_HF.jpg'];   % save
%     saveas(H,fnm)
%     close(H)
% end

% View light-triggered raster and PSTH
TrigEvent = 'BurstOn';
SEvent = 'BurstOff';
FNum = 2;
win = [-0.2 0.5];
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

% Cluster quality
BatchSessionClust(fullpth)

% Behavior
auditory_gonogo_psychplot2(animalID,sessionID)
% auditory_gonogo_psychplot2(animalID,sessionID,[],(1:150))
H = gcf;
fnm = [resdir2 sessionID '_PSYCHPLOT.jpg'];   % save
saveas(H,fnm)
fnm = [resdir2 sessionID '_PSYCHPLOT.fig'];
saveas(H,fnm)

auditory_gonogo_psychplot3(animalID,sessionID)
H = gcf;
maximize_figure(H)
fnm = [resdir2 sessionID '_PSYCHPLOT2.jpg'];   % save
saveas(H,fnm)
fnm = [resdir2 sessionID '_PSYCHPLOT2.fig'];
saveas(H,fnm)