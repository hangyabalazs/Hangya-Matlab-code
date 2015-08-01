function quickanalysis2(animalNO,sessionID,sessionspec)
%QUICKANALYSIS2   Analysis of tetrode data.
%   QUICKANALYSIS2 is designed as an offline analysis tool for tetrode data
%   and behavior, which can be executed in an unsupervised manner on a
%   daily bases. It gives a quick overview of the experiment including
%   response profiles of clustered neurons, light-triggered PSTH and
%   psychometric plots of behavior. It relies on CellBase data handling
%   system.
%
%   QUICKANALYSIS2(ANIMALNO,SESSIONID,SESSIONSPEC) performs the analysis
%   for a session specified by the first two input arguments. SESSIONSPEC
%   should be a 1x3 logical array indicating the need for behavior,
%   recording and stimulation analysis.
%
%   See also INITCB.

% Input argument check
error(nargchk(0,3,nargin))

% Behavior, recording or both
if nargin < 3
    isbeh = 0;
    isrec = 1;
    isstim = 1;
else
    isbeh = sessionspec(1);
    isrec = sessionspec(2);
    isstim = sessionspec(3);
end

% Animal, session
if nargin < 2
    sessionID = '120301b';
end
if nargin < 1
    animalID2 = 'nb029';
    animalID = 'n029';
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

% Position data
% [TimeStamps, ExtractedX, ExtractedY, ExtractedAngle, Targets, Points, NlxHeader] =...
%     Nlx2MatVT(fn,[1 1 1 1 1 1],1,1);

% Convert events file
if isrec
    nlxcsc2mat2(fullpth,'Channels','Events')
end

% Create trial events structure
if isbeh
    TE = solo2trialevents4_auditory_gonogo([fullpth 'data_@auditory_gonogo_balazs_' animalID2 '_' sessionID '.mat']);
    if isrec
        MakeTrialEvents2_gonogo(fullpth)  % synchronize
    end
    
    % Define events
    [ev, ep] = defineEventsEpochs_gonogo;
end

% Update CellBase
if isrec
    addnewsessions
    addnewcells
    cellids = findcell('rat',animalID,'session',sessionID);
    disp(cellids)
end

% Response profiles
if isbeh && isrec
    % Prealign spikes for trial events
    problem_behav_cellid = [];
    for iC = 1:length(cellids),
        cellid = cellids(iC);
        try
            prealignSpikes(cellid,'FUNdefineEventsEpochs',@defineEventsEpochs_gonogo,'filetype','event','ifsave',1,'ifappend',0)
        catch
            disp('Error in prealignSpikes.');
            problem_behav_cellid = [problem_behav_cellid cellid];
        end
    end
    
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
    
%     Lickraster
%     for k = 1:length(cellids)
%         H = figure;
%         viewcell2b(cellids(k),'TriggerName','LickIn','SortEvent','StimulusOff','eventtype','behav','ShowEvents',{{'StimulusOff'}},'Partitions','#ResponseType','window',[-5 5])
%         maximize_figure(H)
%         
%         fnm = [resdir cellid '_HF.jpg'];   % save
%         saveas(H,fnm)
%         close(H)
%     end
end

% Light effects
if isrec && isstim
    
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
end

% Cluster quality
if isrec
    BatchSessionClust(fullpth)
end

% Behavior
if isbeh
    auditory_gonogo_psychplot2(animalID,sessionID)
    % auditory_gonogo_psychplot2(animalID,sessionID,[],(1:150))
    H = gcf;
    fnm = [resdir2 sessionID '_PSYCHPLOT.jpg'];   % save
    saveas(H,fnm)
    fnm = [resdir2 sessionID '_PSYCHPLOT.fig'];
    saveas(H,fnm)
    
    H = auditory_gonogo_psychplot3(animalID,sessionID);
    maximize_figure(H)
    fnm = [resdir2 sessionID '_PSYCHPLOT2.jpg'];   % save
    saveas(H,fnm)
    fnm = [resdir2 sessionID '_PSYCHPLOT2.fig'];
    saveas(H,fnm)
end