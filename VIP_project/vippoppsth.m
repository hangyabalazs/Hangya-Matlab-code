function vippoppsth
%VIPPOPPSTH   Average PSTH.
%   VIPPOPPSTH calculates and plots average adaptive PSTHs for activated,
%   inhibited and tagged group of cells.
%
%   See also VIPPSTH.

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
tagged = find(tagged)';

% Check which Cellbase is active
csb = whichcb;
if ~isequal(csb,'VIP')
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
plot(time,mean(allpsth_act),'Color',[0.32 0.19 0.19])
hold on
errorshade(time,mean(allpsth_act),std(allpsth_act)/sqrt(size(allpsth_act,1)),...
    'LineColor',[0.32 0.19 0.19],'ShadeColor',[0.32 0.19 0.19])
plot(time,mean(allpsth_inh),'Color',[0.48 0.06 0.89])
errorshade(time,mean(allpsth_inh),std(allpsth_inh)/sqrt(size(allpsth_inh,1)),...
    'LineColor',[0.48 0.06 0.89],'ShadeColor',[0.48 0.06 0.89])
plot(time,mean(allpsth_tagged),'Color',[0 0.7 0])
errorshade(time,mean(allpsth_tagged),std(allpsth_tagged)/sqrt(size(allpsth_tagged,1)),...
    'LineColor',[0 0.7 0],'ShadeColor',[0 0.7 0])
line([0 0],[-1.5 6],'Color','k','LineStyle',':')
line([-20 60],[0 0],'Color','k')
set(gca,'box','off',...
    'FontSize',16,'TickDir','out','XLim',[-16 56],'Ylim',[-1.5 6])
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

% Standardize
allpsth2 = allpsth;
for k = 1:size(allpsth,1)
    allpsth2(k,:) = standardize(allpsth(k,:));
end

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

% PSTH
% [psth, spsth, spsth_se] = binraster2psth(spt,dt,sigma,COMPTRIALS,valid_trials);
[psth, spsth, spsth_se] = binraster2apsth(spt,dt,time,sigma,COMPTRIALS,valid_trials);