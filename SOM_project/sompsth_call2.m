function sompsth_call2(I)
%SOMPSTH_CALL2   Creates spike density function with variable Gaussian Kernel.
%   SOMPSTH_CALL2 calls STIMES2BINRASTER for all cells in CellBase (see
%   also CellBase documentation), centered on 'BurstOn' events. Then, spike
%   time series are convolved with a variable width gaussian kernel.
%   Variance of the Gaussian is adapted to the local estimate of spiking
%   probability to implement stronger smoothing when information is sparse.
%
%   SOMPSTH_CALL2(I) runs only for the index set I.
%
%   See also STIMES2BINRASTER.

% Pass the control to the user in case of error
dbstop if error

% Directories
global DATAPATH
fs = filesep;
resdir = [DATAPATH 'SOM' fs 'PSTH' fs];

% Load CellBase
loadcb

% Input argument check
nmc = length(CELLIDLIST);
if nargin < 1
    I = 1:nmc;
end

% Plot raster and PSTH
for k = I;
    disp(k)
    try
        cellid = CELLIDLIST{k};
        [H1 H2] = lpsth(cellid);
        ff1 = [resdir 'SDF_' cellid '.fig'];
        ff2 = [resdir 'PSTH_' cellid '.fig'];
%         saveas(H1,ff1)          % save
%         saveas(H2,ff2)
        close all
    catch ME
        disp(['No raster plot for cell ' num2str(k) '.'])
        disp(ME.message)
    end
end

% -------------------------------------------------------------------------
function [H1 H2] = lpsth(cellid)
% Call STIMES2BINRASTER and calculate adaptive spike density function

% Set input arguments for STIMES2BINRASTER
win = [-0.6 0.6];
dt = 0.001;   % resolution of bin raster in s
EventName = 'PulseOn';
ST = loadcb(cellid,'STIMSPIKES');    % load variables from CellBase
TE = loadcb(cellid,'StimEvents');
epoch_pos1 = findcellstr(ST.events(:,1),EventName);
if epoch_pos1 == 0
    error('Epoch name not found');
end

% Time variables and valid trials
stimes = ST.event_stimes{epoch_pos1};
time = win(1):dt:win(end);
valid_trials = find(~isnan(getfield(TE,EventName)));

% Calculate bin raster
spt = stimes2binraster(stimes(valid_trials),time,dt);

% Trial number and epoch length
[tno tl] = size(spt);

% Calculate adaptive SDF with variable Gaussian Kernel
agvd = zeros(tno,tl);
prob = sum(spt) / tno / (dt * 1000);
figure
A = axes;
hold on
for k = 1:tno   % convolve trial-wise
    spks = find(spt(k,:));
    spno = length(spks);
    for t = 1:spno
        spi = spks(t);
        tspt = zeros(1,tl);
        tspt(spi) = 1;
        if prob(spi) > 1
            keyboard
        end
        wbh = gausswin(9,prob(spi)*50);
        wbh = wbh / sum(wbh);
        plot(wbh,'Color',rand(1,3))
%         ptc = conv(tspt,wbh);
%         lag = (length(ptc)-tl+1) / 2;
%         ptc = ptc(lag:end-lag);
%         agvd(k,:) = agvd(k,:) + ptc;
        agvd(k,:) = agvd(k,:) + filtfilt(wbh,1,tspt);
    end
end
H1 = figure;
subplot(211)
imagesc(agvd);
% set(gca,'XTick',[])
psth_aconv = sum(agvd) / tno / dt;
subplot(212)
plot(time,psth_aconv)
xlim([time(1) time(end)])
axes('Position',[0.7 0.35 0.2 0.2])
Ls = findobj(allchild(A),'type','line');   % copyobj runs to Matlab bug
plot(cell2mat(get(Ls,'XData'))',cell2mat(get(Ls,'YData'))')

% Input variables for 'viewcell2b
SEvent = 'BurstOff';
FNum = 2;
parts = 'all';
sigma = 0.001;
PSTHstd = 'on';
ShEvent = {{'BurstOff'}};
ShEvColors = hsv(length(ShEvent{1}));
ShEvColors = mat2cell(ShEvColors,ones(size(ShEvColors,1),1),3);

% Call 'viewcell2b' for comparison
H2 = figure;
set(gcf,'renderer','painters')   % temporaray change renderer because OpenGL locks the plot which result an error in legend layout handling
viewcell2b(cellid,'TriggerName',EventName,'SortEvent',SEvent,'ShowEvents',ShEvent,'ShowEventsColors',{ShEvColors},...
    'FigureNum',FNum,'eventtype','stim','window',win,'dt',dt,'sigma',sigma,'PSTHstd',PSTHstd,'Partitions',parts,...
    'EventMarkerWidth',0,'PlotZeroLine','on')
pause(0.05)   % if reset the renderer two early, the same error occurs
set(gcf,'renderer','opengl')   % reset renderer