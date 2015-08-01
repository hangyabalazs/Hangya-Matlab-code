function vippoppsth_A1
%VIPPOPPSTH_A1   Average PSTH.
%   VIPPOPPSTH_A1 calculates and plots average adaptive PSTHs for activated,
%   inhibited and tagged group of cells.
%
%   See also VIPPSTH.

% Load
load('C:\Balazs\_analysis\VIP\A1_psth_variables.mat')

% Groups
frlim = 1;   % lower firing rate limit for detecting inhibition
tagged = logical(vldty) & (isact==2);   % tagged cells
inx_act = logical(vldty) & (isact==1) & (isact~=2);  % activated cells
inx_inh = logical(vldty) & (isinh) & (baseline>frlim) & (isact~=2);   % inhibited cells; firing rate > 1Hz criterion
activated = find(inx_act&~inx_inh);  % indices of activated only cells
inhibited = find(inx_inh&~inx_act);  % indices of inhibited only cells
ai = find(inx_act&inx_inh);   % indices of cells with dual effect
inx = activation_peak(ai) < inhibition_peak(ai);
activated_inhibited = sort(ai(inx));  % activated, then inhibited
inhibited_activated = ai(~inx);   % inhibited, then activated
activated = [activated; activated_inhibited]';   % activated and activated-inhibited
inhibited = [inhibited; inhibited_activated]';   % inhibited and inhibited-activated
tagged = find(tagged)';   % indices of tagged cells

% Check which Cellbase is active
csb = whichcb;
if ~isequal(csb,'VIP_A1')
    disp('Switch to VIP CellBase!')
    return
end

% PSTH
[allpsth_act wn] = poppsth(activated);
[allpsth_inh wn] = poppsth(inhibited);
[allpsth_tagged wn] = poppsth(tagged);
time = linspace(wn(1)*1000,wn(2)*1000,size(allpsth_act,2));

% Plot
figure
plot(time,mean(allpsth_act),'Color',[1 0.55 0.33])
hold on
errorshade(time,mean(allpsth_act),std(allpsth_act)/sqrt(size(allpsth_act,1)),...
    'LineColor',[1 0.55 0.33],'ShadeColor',[1 0.55 0.33])
plot(time,mean(allpsth_inh),'Color',[0.48 0.06 0.89])
errorshade(time,mean(allpsth_inh),std(allpsth_inh)/sqrt(size(allpsth_inh,1)),...
    'LineColor',[0.48 0.06 0.89],'ShadeColor',[0.48 0.06 0.89])
plot(time,mean(allpsth_tagged),'Color',[0 0.7 0])
errorshade(time,mean(allpsth_tagged),std(allpsth_tagged)/sqrt(size(allpsth_tagged,1)),...
    'LineColor',[0 0.7 0],'ShadeColor',[0 0.7 0])
line([0 0],[-4 8],'Color','k','LineStyle',':')
line([-200 600],[0 0],'Color','k')
set(gca,'box','off',...
    'FontSize',16,'TickDir','out','XLim',[-200 600],'Ylim',[-2.6 8])
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
    dt = 0.002;   % resolution of bin raster in s (was 5 ms before for A1)
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