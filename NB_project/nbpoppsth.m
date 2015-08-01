function nbpoppsth(I,issave)
%NBPOPPSTH   Average PSTH.
%   NBPOPPSTH(I,ISSAVE) calculates and plots average adaptive PSTHs for a
%   set of cells (see APSTH). Input parameters: 
%       I - index set to CELLIDLIST (see CellBase documentation); if empty
%           or not specified, all well-separated cells are selected (ID>20,
%           L-ratio<0.15; see LRATIO)
%       ISSAVE - controls saving
%
%   See also APSTH, LRATIO and NBPOPPSTH_CALL.

%   Edit log: BH 7/3/12
 
% Pass the control to the user in case of error
dbstop if error
 
% Input argument check
error(nargchk(0,2,nargin))
if nargin < 2
    issave = false;
end
if nargin < 1
    I = [];
end
 
% Directories
global DATAPATH
fs = filesep;
resdir = [DATAPATH 'NB' fs 'poppsth' fs];
 
% List of cellIDs
if isempty(I)
    Lratio = getvalue('Lr_PC');
    ID = getvalue('ID_PC');
    ptinx = ID > 20 & Lratio < 0.15;
    I = find(ptinx);
end

% PSTH
[allpsth_orig allpsth allspsth_orig allspsth wn] = poppsth(I);
time = linspace(wn(1)*1000,wn(2)*1000,size(allpsth,2));

% Plot
figure
plot(time,mean(allpsth),'Color',[0.32 0.19 0.19])
hold on
errorshade(time,mean(allpsth),std(allpsth)/sqrt(size(allpsth,1)),...
    'LineColor',[0.32 0.19 0.19],'ShadeColor',[0.32 0.19 0.19])
% line([0 0],[-1.5 6],'Color','k','LineStyle',':')
% line([-20 60],[0 0],'Color','k')
% set(gca,'box','off',...
%     'FontSize',16,'TickDir','out','XLim',[-16 56],'Ylim',[-1.5 6])
xlabel('Time')
ylabel('Normalized firing rate')
keyboard

% [mx mxinx] = max(allpsth,[],2);
% [srt srtinx] = sort(mxinx,'ascend');
% figure
% imagesc(time,1:19,allpsth(srtinx,:))
figure
imagesc(time,1:19,allpsth)
xlim([-600 600])
colormap(hot)
keyboard

% -------------------------------------------------------------------------
function [allpsth allpsth2 allspsth allspsth2 wn] = poppsth(I)

% Load CellBase
loadcb

% Call 'main'
allpsth = [];
allspsth = [];
for k = I
    cellid = CELLIDLIST{k}; %#ok<USENS>
    disp(cellid)
    [psth spsth wn] = main(cellid);
    allpsth = [allpsth; psth]; %#ok<AGROW>
    allspsth = [allspsth; spsth]; %#ok<AGROW>
end

% Standardize
allpsth2 = allpsth;
allspsth2 = allspsth;
for k = 1:size(allpsth,1)
    allpsth2(k,:) = standardize(allpsth(k,:));
    allspsth2(k,:) = standardize(allspsth(k,:));
end

% -------------------------------------------------------------------------
function [psth spsth win] = main(cellid)

% Input argument check
if nargin < 2
    win = [-0.65 0.65];  % time window for bin raster
    dt = 0.004;   % resolution of bin raster in s
    dsply = 0;   % supress display
    sigma = 0.02;
    parts = 'all';
end

% Set parameters and load CellBase variables
EventName = 'LeftPortIn';
ST = loadcb(cellid,'EVENTSPIKES');   % load prealigned spikes for stimulation events
TE = loadcb(cellid,'TrialEvents');
event_pos = findcellstr(ST.events(:,1),EventName);
if event_pos == 0
    error('Event name not found');
end
stimes = ST.event_stimes{event_pos};
time = win(1):dt:win(end);
valid_trials = find(~isnan(TE.Hit));    % Hits only

% Calculate bin rasters
% spt = stimes2binraster(stimes1(valid_trials),time,dt);
spt = stimes2binraster(stimes,time,dt);

% Partition trials
[COMPTRIALS, TAGS] = partition_trials(TE,parts);

% PSTH
% [psth, spsth, spsth_se] = binraster2psth(spt,dt,sigma,COMPTRIALS,valid_trials);
[psth, spsth, spsth_se] = binraster2apsth(spt,dt,time,sigma,COMPTRIALS,valid_trials);