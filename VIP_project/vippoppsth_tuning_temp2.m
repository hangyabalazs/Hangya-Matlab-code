function vippoppsth_tuning
%VIPPOPPSTH_TUNING   Average PSTH.
%   VIPPOPPSTH_TUNING calculates and plots average adaptive PSTHs.
%
%   See also VIPPSTH and VIPPOPPSTH2.

% Check which Cellbase is active
csb = whichcb;
if ~isequal(csb,'VIP_A1')
    disp('Switch to VIP CellBase!')
    return
end

% Groups
Sound = [28 46 56 99 126 146 162 167 168 169 208 212 216 218 223 234 238 ...
    250 259 271 279 282 288 303 305 311 314 319 320 329 330 353 357 363];
NoSound = [49 55 57 60 63 65 66 70 149 171 173 177 276 286 299 346 359];

% PSTH
[allpsth_s wn] = poppsth(Sound);
[allpsth_ns wn] = poppsth(NoSound);
time = linspace(wn(1)*1000,wn(2)*1000,size(allpsth_s,2));

% Plot
figure
plot(time,mean(allpsth_s),'Color',[0.6 0.2 0])
hold on
errorshade(time,mean(allpsth_s),std(allpsth_s)/sqrt(size(allpsth_s,1)),...
    'LineColor',[0.6 0.2 0],'ShadeColor',[0.6 0.2 0])
plot(time,mean(allpsth_ns),'Color',[0.87 0.47 0])
errorshade(time,mean(allpsth_ns),std(allpsth_ns)/sqrt(size(allpsth_ns,1)),...
    'LineColor',[0.87 0.47 0],'ShadeColor',[0.87 0.47 0])
line([0 0],[-1.5 6],'Color','k','LineStyle',':')
line([-200 600],[0 0],'Color','k')
set(gca,'box','off',...
    'FontSize',16,'TickDir','out','XLim',[-200 600],'Ylim',[-1 4.5])
xlabel('Time')
ylabel('Normalized firing rate')
keyboard

% -------------------------------------------------------------------------
function [allpsth2 wn] = poppsth(I)

% Pass the control to the user in case of error
dbstop if error

% Directories
global DATAPATH
fs = filesep;
resdir = [DATAPATH 'VIP' fs 'poppsth_A1_' fs];

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

% Standardize
allpsth2 = allpsth;
for k = 1:size(allpsth,1)
    allpsth2(k,:) = standardize(allpsth(k,:));
end

% -------------------------------------------------------------------------
function [spsth win] = main(cellid)

% Input argument check
if nargin < 2
    win = [-0.25 0.65];  % time window for bin raster
    dt = 0.005;   % resolution of bin raster in s
    dsply = 0;   % supress display
    sigma = 0.001;
    parts = 'all';
end

% Set parameters and load CellBase variables
EventName1 = 'PulseOn';
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

% PSTH
% [psth, spsth, spsth_se] = binraster2psth(spt,dt,sigma,COMPTRIALS,valid_trials);
[psth, spsth, spsth_se] = binraster2apsth(spt,dt,sigma,COMPTRIALS,valid_trials);