function sompsth_call(I)
%SOMPSTH_CALL   Creates raster plot and PSTH.
%   SOMPSTH_CALL calls VIEWCELL2B for all cells in CellBase (see also
%   CellBase documentation). A plot whith spike raster and peri-stimulus
%   time histogram is created, centered on 'BurstOn' events.
%
%   SOMPSTH_CALL(I) runs only for the index set I.
%
%   See also VIEWCELL2B.

% Pass the control to the user in case of error
dbstop if error

% Directories
global DATAPATH
fs = filesep;
resdir = [DATAPATH 'SOM' fs 'KL' fs];

% Load CellBase
loadcb

% Input argument check
nmc = length(CELLIDLIST);
if nargin < 1
    I = 1:nmc;
end

% Plot raster and PSTH
ef = nan(1,nmc);
ff1 = [resdir 'effect.mat'];
for k = I;
    disp(k)
    try
        figure
        lpsth(CELLIDLIST{k})
        ef(k) = input('Inhibited (-1), excited (1), neither (0) or both(2)? ');
        close all
%         save(ff1,'ef')          % save
    catch ME
        disp(['No raster plot for cell ' num2str(k) '.'])
        disp(ME.message)
    end
end

% -------------------------------------------------------------------------
function lpsth(cellid)
% Call VIEWCELL2B to get raster plot and PSTH

% Set input arguments for VIEWCELL2B
win = [-0.6 0.6];
dt = 0.001;   % resolution of bin raster in s
EventName = 'BurstOn';
ST = loadcb(cellid,'STIMSPIKES');    % load variables from CellBase
TE = loadcb(cellid,'StimEvents');

epoch_pos1 = findcellstr(ST.events(:,1),EventName);
if epoch_pos1 == 0
    error('Epoch name not found');
end

SEvent = 'BurstOff';
FNum = 2;
parts = 'all';
sigma = 0.001;
PSTHstd = 'on';
ShEvent = {{'BurstOff'}};
ShEvColors = hsv(length(ShEvent{1}));
ShEvColors = mat2cell(ShEvColors,ones(size(ShEvColors,1),1),3);

% Call 'viewcell2b'
set(gcf,'renderer','painters')   % temporaray change renderer because OpenGL locks the plot which result an error in legend layout handling
viewcell2b(cellid,'TriggerName',EventName,'SortEvent',SEvent,'ShowEvents',ShEvent,'ShowEventsColors',{ShEvColors},...
    'FigureNum',FNum,'eventtype','stim','window',win,'dt',dt,'sigma',sigma,'PSTHstd',PSTHstd,'Partitions',parts,...
    'EventMarkerWidth',0,'PlotZeroLine','on')
pause(0.05)   % if reset the renderer two early, the same error occurs
set(gcf,'renderer','opengl')   % reset renderer