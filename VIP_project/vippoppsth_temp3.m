function vippoppsth

% Load
load('C:\Balazs\_analysis\VIP\psth_variables.mat')

% Groups
tagged =  logical(vldty) & (isact==2);
inx_act = logical(vldty) & (isact==1) & (baseline>2);
inx_inh = logical(vldty) & (isinh) & (baseline>2);
activated = setdiff(find(inx_act&~inx_inh),tagged);
inhibited = find(inx_inh&~inx_act);
ai = find(inx_act&inx_inh);
inx = activation_peak(ai) < inhibition_peak(ai);
activated_inhibited = sort([ai(inx); 10]);    % 10 is activated-inhibited
inhibited_activated = setdiff(ai(~inx),10);
activated = [activated; activated_inhibited]';
inhibited = [inhibited; inhibited_activated]';

% Check which Cellbase is active
csb = whichcb;
if ~isequal(csb,'VIP')
    disp('Switch to VIP CellBase!')
    return
end

allpsth_act = poppsth(activated);
allpsth_inh = poppsth(inhibited);


figure
plot(linspace(wn(1)*1000,wn(2)*1000,size(allpsth2,2)),mean(allpsth2))
xlim([-20 55])
% line([0 0],[-1.5 1],'Color','k','LineWidth',2)
line([0 0],[-1 2],'Color','k','LineWidth',2)
box off

% -------------------------------------------------------------------------
function poppsth(I)

% Pass the control to the user in case of error
dbstop if error

% Directories
global DATAPATH
fs = filesep;
resdir = [DATAPATH 'VIP' fs 'poppsth' fs];

% Load CellBase
loadcb

% Input argument check
nmc = length(CELLIDLIST);
if nargin < 1
    I = 1:nmc;
end

% Call 'main'
allpsth = [];
for k = I
    disp(k)
    cellid = CELLIDLIST{k};
%     try
        [psth wn] = main(cellid);
        allpsth = [allpsth; psth];
%     catch ME
%         disp(['No PSTH calculated for cell ' num2str(k) '.'])
%         disp(ME.message)
%     end
end
keyboard

% Standardize
allpsth2 = allpsth;
for k=1:size(allpsth,1)
    allpsth2(k,:) = standardize(allpsth(k,:));
end

% Plot PSTH for activated ones
uiopen('C:\Balazs\_analysis\VIP\vippsth\PSTH_h016_101128a_6.3.fig',1)       % inh?
uiopen('C:\Balazs\_analysis\VIP\vippsth\PSTH_h016_101203a_2.1.fig',1)
uiopen('C:\Balazs\_analysis\VIP\vippsth\PSTH_h023_110111a_3.1.fig',1)
uiopen('C:\Balazs\_analysis\VIP\vippsth\PSTH_h023_110114a_5.2.fig',1)
uiopen('C:\Balazs\_analysis\VIP\vippsth\PSTH_h023_110118a_2.1.fig',1)       % inh?
uiopen('C:\Balazs\_analysis\VIP\vippsth\PSTH_h026_110221a_1.1.fig',1)
uiopen('C:\Balazs\_analysis\VIP\vippsth\PSTH_h058_110729a_6.1.fig',1)
uiopen('C:\Balazs\_analysis\VIP\vippsth\PSTH_h058_110802a_3.1.fig',1)
uiopen('C:\Balazs\_analysis\VIP\vippsth\PSTH_h016_101203b_6.1.fig',1)
uiopen('C:\Balazs\_analysis\VIP\vippsth\PSTH_h023_110116a_5.2.fig',1)
uiopen('C:\Balazs\_analysis\VIP\vippsth\PSTH_h058_110804a_3.1.fig',1)
uiopen('C:\Balazs\_analysis\VIP\vippsth\PSTH_h058_110810a_1.1.fig',1)

% Plot PSTH for inhibited ones
uiopen('C:\Balazs\_analysis\VIP\vippsth\PSTH_h023_110122a_3.1.fig',1)
uiopen('C:\Balazs\_analysis\VIP\vippsth\PSTH_h026_110220a_1.2.fig',1)
uiopen('C:\Balazs\_analysis\VIP\vippsth\PSTH_h058_110603a_4.7.fig',1)
uiopen('C:\Balazs\_analysis\VIP\vippsth\PSTH_h058_110608a_4.1.fig',1)
uiopen('C:\Balazs\_analysis\VIP\vippsth\PSTH_h058_110731a_4.2.fig',1)
uiopen('C:\Balazs\_analysis\VIP\vippsth\PSTH_h058_110802a_5.2.fig',1)
uiopen('C:\Balazs\_analysis\VIP\vippsth\PSTH_h016_101203b_1.1.fig',1)       % act?
uiopen('C:\Balazs\_analysis\VIP\vippsth\PSTH_h058_110603a_6.1.fig',1)
uiopen('C:\Balazs\_analysis\VIP\vippsth\PSTH_h058_110606a_4.1.fig',1)
uiopen('C:\Balazs\_analysis\VIP\vippsth\PSTH_h058_110606a_6.1.fig',1)
uiopen('C:\Balazs\_analysis\VIP\vippsth\PSTH_h058_110609a_6.3.fig',1)
uiopen('C:\Balazs\_analysis\VIP\vippsth\PSTH_h058_110731a_5.1.fig',1)
uiopen('C:\Balazs\_analysis\VIP\vippsth\PSTH_h058_110731a_6.1.fig',1)
uiopen('C:\Balazs\_analysis\VIP\vippsth\PSTH_h058_110802a_5.1.fig',1)
uiopen('C:\Balazs\_analysis\VIP\vippsth\PSTH_h058_110807a_5.1.fig',1)       % act?

uiopen('C:\Balazs\_analysis\VIP\vippsth\PSTH_h016_101203b_1.1.fig',1)
uiopen('C:\Balazs\_analysis\VIP\vipisinfluenced11\SDF_h016_101203b_1.1.fig',1)
uiopen('C:\Balazs\_analysis\VIP\vipisinfluenced11\RASTER_h016_101203b_1.1.fig',1)

% -------------------------------------------------------------------------
function [spsth win] = main(cellid)

% Input argument check
if nargin < 2
    win = [-0.02 0.06];  % time window for bin raster
    dt = 0.001;   % resolution of bin raster in s
    dsply = 0;   % supress display
    sigma = 0.001;
    parts = 'all';
end

% Set parameters and load CellBase variables
EventName1 = 'BurstOn';
ST = loadcb(cellid,'STIMSPIKES');   % load prealigned spikes for stimulation events
TE = loadcb(cellid,'StimEvents');
epoch_pos1 = findcellstr(ST.events(:,1),EventName1);
if epoch_pos1 == 0
    error('Epoch name not found');
end
stimes1 = ST.event_stimes{epoch_pos1};
time = win(1):dt:win(end);
valid_trials = find(~isnan(getfield(TE,EventName1)));

% Calculate bin rasters
% spt = stimes2binraster(stimes1(valid_trials),time,dt);
spt = stimes2binraster(stimes1,time,dt);

% Partition trials
[COMPTRIALS, TAGS] = partition_trials(TE,parts);

[psth, spsth, spsth_se] = binraster2psth(spt,dt,sigma,COMPTRIALS,valid_trials);
[psth, spsth, spsth_se] = binraster2apsth(spt,dt,time,sigma,COMPTRIALS,valid_trials);