%% stim raster

% View light-triggered raster and PSTH
TrigEvent = 'PulseOn';
SEvent = 'BurstOff';
win = [-0.2 0.5];
parts = 'all';
% parts = '#BurstNPulse';
dt = 0.001;
sigma = 0.001;
PSTHstd = 'on';
ShEvent = {{'PulseOn'}};
ShEvColors = hsv(length(ShEvent{1}));
ShEvColors = mat2cell(ShEvColors,ones(size(ShEvColors,1),1),3);
H = figure;
viewcell2b(cellid,'TriggerName',TrigEvent,'SortEvent',SEvent,'ShowEvents',ShEvent,'ShowEventsColors',{ShEvColors},...
    'FigureNum',H,'eventtype','stim','window',win,'dt',dt,'sigma',sigma,'PSTHstd',PSTHstd,'Partitions',parts,...
    'EventMarkerWidth',0,'PlotZeroLine','off')

%%

wn = [-0.01 0.02];
dt = 0.001;
[psth, spsth, spsth_se, ~, spt] = ...
    ultimate_psth(cellid,'stim','BurstOn',wn,...
    'dt',dt,'display',true,'sigma',0.001,'parts','all','isadaptive',1,...
    'event_filter','custom','filterinput','BurstNPulse==20','maxtrialno',Inf);
sspt = sum(spt);
figure
bar(wn(1):dt:wn(2),sspt,'FaceColor','k','EdgeColor','k','BarWidth',1)
xlim(wn)

%% hit raster

H = figure;
pause(0.01)
viewcell2b(cellid,'TriggerName','LeftPortIn','SortEvent',...
    'StimulusOn','eventtype','behav','ShowEvents',{{'StimulusOn'}},...
    'Partitions','#ResponseType','window',[-0.5 1])

%% auto-correlation

acg(cellid)

%% cross-correlation

cellid1 = 'n038_120922a_3.5';
cellid2 = 'n038_120922a_5.1';
ccg_Tourette({cellid1 cellid2})