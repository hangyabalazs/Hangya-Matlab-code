function viprebound
%VIPREBOUND   Average PSTH aligned to inhibition offset.
%   VIPREBOUND calculates and plots average adaptive PSTHs for the 
%   inhibited-activated group of cells aligned to the first spike after the 
%   trough of inhibition.
%
%   See also VIPPOPPSTH_A1.

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
[allpsth_inhact wn] = poppsth(inhibited_activated,inhibition_peak(inhibited_activated));
time = linspace(wn(1)*1000,wn(2)*1000,size(allpsth_inhact,2));

% Plot
figure
hold on
errorshade(time,mean(allpsth_inhact),std(allpsth_inhact)/sqrt(size(allpsth_inhact,1)),...
    'LineColor',[0.48 0.06 0.89],'ShadeColor',[0.48 0.06 0.89])
line([0 0],[-4 8],'Color','k','LineStyle',':')
line([-200 600],[0 0],'Color','k')
set(gca,'box','off',...
    'FontSize',16,'TickDir','out','XLim',[-200 600],'Ylim',[-2.6 8])
xlabel('Time')
ylabel('Normalized firing rate')
keyboard

% -------------------------------------------------------------------------
function [allpsth2 wn] = poppsth(I,inhibition_peak)

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
I = I(:)';   % convert to row vector

% Call 'main'
allpsth = [];
for k = 1:length(I)
    disp(I(k))
    cellid = CELLIDLIST{I(k)};
%     try
        [psth wn] = main(cellid,inhibition_peak(k)/1000);    % inhibition peak time is converted to s from ms
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
function [spsth win] = main(cellid,inhibition_peak)

% Input argument check
if nargin < 3
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

% Realign spikes
NumTrials = length(stimes1);   % number of trials
for iTrial = 1:NumTrials   % loop through all trials
    cst = stimes1{iTrial};   % current trial
    cst = cst - inhibition_peak;  % align to inhibition peak
    fsp = min(cst(cst>0));   % first spike after inhibition_peak
    if(~isempty(fsp))   % if there was a spike after the peak of the inhibition
        cst = cst - fsp;   % align to first spike after inhibition_peak
        cst(cst==0) = [];   % omit the spike to which we aligned
    end
    stimes1{iTrial} = cst;
end

% Calculate bin rasters
% spt = stimes2binraster(stimes1(valid_trials),time,dt);
spt = stimes2binraster(stimes1,time,dt);

% Partition trials
[COMPTRIALS, TAGS] = partition_trials(TE,parts);

% PSTH
% [psth, spsth, spsth_se] = binraster2psth(spt,dt,sigma,COMPTRIALS,valid_trials);
[psth, spsth, spsth_se] = binraster2apsth(spt,dt,sigma,COMPTRIALS,valid_trials);