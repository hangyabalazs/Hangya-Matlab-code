function quickanalysis2_replot(animalNO,sessionID,sessionspec,protocoltag)
%QUICKANALYSIS2_REPLOT   Analysis of tetrode data.
%   QUICKANALYSIS2_REPLOT is designed as an offline analysis tool for
%   tetrode data and behavior, which can be executed in an unsupervised
%   manner on a daily bases. It gives a quick overview of the experiment
%   including response profiles of clustered neurons, light-triggered PSTH
%   and psychometric plots of behavior. It relies on CellBase data handling
%   system.
%
%   QUICKANALYSIS2_REPLOT(ANIMALNO,SESSIONID,SESSIONSPEC) performs the
%   analysis for a session specified by the first two input arguments.
%   SESSIONSPEC should be a 1x3 logical array indicating the need for
%   behavior, recording and stimulation analysis.
%
%   QUICKANALYSIS2_REPLOT(ANIMALNO,SESSIONID,SESSIONSPEC,PROTOCOLTAG)
%   accepts a PROTOCOLTAG argument to allow calls to trial event conversion
%   programs for different versions of the behavioral protocol.
%
%   QUICKANALYSIS2_REPLOT does not do file conversion and spike
%   prealignment, assuming those steps had already been performed. It plots
%   and saves PSTHs.
%
%   See also QUICKANALYSIS2.

% Input argument check
error(nargchk(0,4,nargin))

% Behavior, recording or both
if nargin < 4
    protocoltag = '';
end
if nargin < 3
    isbeh = 1;
    isrec = 1;
    isstim = 1;
else
    isbeh = sessionspec(1);
    isrec = sessionspec(2);
    isstim = sessionspec(3);
end

% Animal, session
if nargin < 2
    sessionID = '121217x';
end
if nargin < 1
    animalID2 = 'nb045';
    animalID = 'n045';
else
    animalID2 = ['nb0' num2str(animalNO)];
    animalID = ['n0' num2str(animalNO)];
end

fullpth = [getpref('cellbase','datapath') '\' animalID '\' sessionID '\'];

% Stop if error
dbstop if error

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

% Find cells
if isrec
    cellids = findcell('rat',animalID,'session',sessionID);
    disp(cellids)
end

% % Response profiles
% if isbeh && isrec
%     
%     % Is predictive?
%     for k = 1:length(cellids)
%         H = figure;
%         viewcell2b(cellids(k),'TriggerName','StimulusOn','SortEvent','TrialStart','eventtype','behav','ShowEvents',{{'StimulusOff'}},'Partitions','#ResponseType','window',[-5 5])
%         maximize_figure(H)
%         
%         cellidt = cellids{k};
%         cellidt(cellidt=='.') = '_';
%         fnm = [resdir cellidt '_IPD.jpg'];   % save
%         saveas(H,fnm)
%         close(H)
%     end
%     
%     % Hit & FA
%     for k = 1:length(cellids)
%         H = figure;
%         pause(0.01)
%         viewcell2b(cellids(k),'TriggerName','LeftPortIn','SortEvent','StimulusOn','eventtype','behav','ShowEvents',{{'StimulusOn'}},'Partitions','#ResponseType','window',[-5 5])
%         maximize_figure(H)
%         
%         cellidt = cellids{k};
%         cellidt(cellidt=='.') = '_';
%         fnm = [resdir cellidt '_HF.jpg'];   % save
%         saveas(H,fnm)
%         close(H)
%     end
%     
% %     % Hit & FA #2
% %     for k = 1:length(cellids)
% %         H = figure;
% %         pause(0.01)
% %         viewcell2b(cellids(k),'TriggerName','DeliverFeedback','SortEvent','PseudoStimulusOn','eventtype','behav','ShowEvents',{{'LeftPortIn'}},'Partitions','#ResponseType','window',[-5 5])
% %         maximize_figure(H)
% %         
% %         cellidt = cellids{k};
% %         cellidt(cellidt=='.') = '_';
% %         fnm = [resdir cellidt '_HF2.jpg'];   % save
% %         saveas(H,fnm)
% %         close(H)
% %     end
%     
%     % Does it depend on stim intensity?
%     for k = 1:length(cellids)
%         H = figure;
%         viewcell2b(cellids(k),'TriggerName','StimulusOn','SortEvent','StimulusOff','eventtype','behav','ShowEvents',{{'StimulusOff'}},'Partitions','#StimulusDuration&Hit','window',[-5 5])
%         maximize_figure(H)
%         
%         cellidt = cellids{k};
%         cellidt(cellidt=='.') = '_';
%         fnm = [resdir cellidt '_SI.jpg'];   % save
%         saveas(H,fnm)
%         close(H)
%     end
% end
%     
% % Light effects
% if isrec && isstim
%     
%     % View light-triggered raster and PSTH
%     TrigEvent = 'BurstOn';
%     SEvent = 'BurstOff';
%     win = [-0.2 0.5];
%     % parts = 'all';
%     parts = '#BurstNPulse';
%     dt = 0.001;
%     sigma = 0.001;
%     PSTHstd = 'on';
%     ShEvent = {{'PulseOn','PulseOff','BurstOff'}};
%     ShEvColors = hsv(length(ShEvent{1}));
%     ShEvColors = mat2cell(ShEvColors,ones(size(ShEvColors,1),1),3);
%     for iCell = 1:length(cellids)
%         cellid = cellids(iCell);
%         H = figure;
%         viewcell2b(cellid,'TriggerName',TrigEvent,'SortEvent',SEvent,'ShowEvents',ShEvent,'ShowEventsColors',{ShEvColors},...
%             'FigureNum',H,'eventtype','stim','window',win,'dt',dt,'sigma',sigma,'PSTHstd',PSTHstd,'Partitions',parts,...
%             'EventMarkerWidth',0,'PlotZeroLine','off')
%         maximize_figure(H)
%         
%         cellidt = cellid{1};
%         cellidt(cellidt=='.') = '_';
%         fnm = [resdir cellidt '_LS.jpg'];   % save
%         saveas(H,fnm)
%         close(H)
%     end
% end
% 
% Cluster quality
if isrec
%     fdir = fullfile(fullpth,'FD');
%     if ~isdir(fdir)
%         mkdir(fdir)
%     end
    BatchSessionClust(fullpth)
end
% 
% % Behavior
% if isbeh
%     auditory_gonogo_psychplot2(animalID,sessionID)
%     % auditory_gonogo_psychplot2(animalID,sessionID,[],(1:150))
%     H = gcf;
%     fnm = [resdir2 sessionID '_PSYCHPLOT.jpg'];   % save
%     saveas(H,fnm)
%     fnm = [resdir2 sessionID '_PSYCHPLOT.fig'];
%     saveas(H,fnm)
%     
%     H = auditory_gonogo_psychplot3(animalID,sessionID);
%     maximize_figure(H)
%     fnm = [resdir2 sessionID '_PSYCHPLOT2.jpg'];   % save
%     saveas(H,fnm)
%     fnm = [resdir2 sessionID '_PSYCHPLOT2.fig'];
%     saveas(H,fnm)
% end