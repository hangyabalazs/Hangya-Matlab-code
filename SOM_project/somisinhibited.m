function p_value = somisinhibited(cellid)
%SOMISINHIBITED   Test for event-locked inhibition.
%   P = SOMISINHIBITED(CELLID) tests whether spike number after an event
%   (e.g. 'PulseOn') decreases. It uses 10 ms time bins for spike number
%   calculation and a before-event baseline for null-hypothesis estimataion
%   of spike number. Instead of taking the first 10 ms after the event, it
%   uses the first after-event spike for starting point of the window given
%   a spike within 10 ms occurs.
%
%   SOMISINHIBITED takes a different bin raster (typically aligned to
%   'BurstOn' events) for baseline distribution than the the one for
%   testing against baseline (typically aligned to 'PulseOn' events).
%
%   See also SOMISSTIM3 and STIMES2BINRASTER.

% Input argument check
if nargin < 2
    win = [-1.0 0.6];  % time window for bin raster
    dt = 0.001;   % resolution of bin raster in s
    dsply = 0;   % repress display
end

% Set parameters and load CellBase variables
EventName1 = 'BurstOn';
EventName2 = 'PulseOn';
ST = loadcb(cellid,'STIMSPIKES');   % load prealigned spikes for stimulation events
TE = loadcb(cellid,'StimEvents');
epoch_pos1 = findcellstr(ST.events(:,1),EventName1);
epoch_pos2 = findcellstr(ST.events(:,1),EventName2);
if epoch_pos1 == 0 || epoch_pos2 == 0
    error('Epoch name not found');
end
stimes1 = ST.event_stimes{epoch_pos1};
stimes2 = ST.event_stimes{epoch_pos2};
time = win(1):dt:win(end);
valid_trials1 = find(~isnan(getfield(TE,EventName1)));
minfreq = min([TE.BurstNPulse]);
maxpow = max([TE.PulsePower]);
inx = ~isnan(getfield(TE,EventName2)) & TE.BurstNPulse==minfreq & TE.PulsePower==maxpow;
valid_trials2 = find(inx);

% Calculate bin rasters
spt1 = stimes2binraster(stimes1(valid_trials1),time,dt);
spt2 = stimes2binraster(stimes2(valid_trials2),time,dt);

% Set input arguments for rater plot and PSTH
if dsply
    SEvent = 'BurstOff';
    FNum = 2;
    parts = 'all';
    sigma = 0.001;
    PSTHstd = 'on';
    ShEvent = {{'BurstOff'}};
    ShEvColors = hsv(length(ShEvent{1}));
    ShEvColors = mat2cell(ShEvColors,ones(size(ShEvColors,1),1),3);
    
    % Plot raster plot and PSTH for 'BurstOn'
    set(gcf,'renderer','painters')   % temporaray change renderer because OpenGL locks the plot which result an error in legend layout handling
    viewcell2b(cellid,'TriggerName',EventName1,'SortEvent',SEvent,'ShowEvents',ShEvent,'ShowEventsColors',{ShEvColors},...
        'FigureNum',FNum,'eventtype','stim','window',win,'dt',dt,'sigma',sigma,'PSTHstd',PSTHstd,'Partitions',parts,...
        'EventMarkerWidth',0,'PlotZeroLine','on')
    pause(0.05)   % if reset the renderer two early, the same error occurs
    set(gcf,'renderer','opengl')   % reset renderer
    
    % Plot raster plot and PSTH for 'PulseOn'
    set(gcf,'renderer','painters')   % temporaray change renderer because OpenGL locks the plot which result an error in legend layout handling
    viewcell2b(cellid,'TriggerName',EventName2,'SortEvent',SEvent,'ShowEvents',ShEvent,'ShowEventsColors',{ShEvColors},...
        'FigureNum',FNum,'eventtype','stim','window',win,'dt',dt,'sigma',sigma,'PSTHstd',PSTHstd,'Partitions',parts,...
        'EventMarkerWidth',0,'PlotZeroLine','on')
    pause(0.05)   % if reset the renderer two early, the same error occurs
    set(gcf,'renderer','opengl')   % reset renderer
end

% Calculate spike numbers and p values
res = 10;   % resolution in ms
dtt = dt * 1000;   % resolution of bin raster in ms
wn = win * 1000;   % window boundaries in ms
p_value = bspn(spt1,spt2,dtt,wn,res);

% -------------------------------------------------------------------------
function p_value = bspn(spt_baseline,spt_test,dt,win,res)

% Trial number and epoch length
[tno tl] = size(spt_baseline);

% Number of bins
nmbn = round(res/dt);

% Pre-stimulus time window to consider for null hypothesis
st = abs(win(1)) / dt;   % number of pre-stim values in 'spt'

% Spike numbers - baseline
nm = floor(st/nmbn);
pspno = nan(tno,nm);    % spike number
next = 1;
for t = 1:nmbn:st
    for k = 1:tno
        cspt = spt_baseline(k,t:t+nmbn-1);
        pspno(k,next) = sum(cspt);
    end
    next = next + 1;
end
spno2 = sum(reshape(spt_baseline(1:tno*nmbn*round(tl/nmbn)),tno*nmbn,round(tl/nmbn)));

% Spike numbers - test
tno_test = size(spt_test,1);
pspno_tt = nan(tno_test,1);
for k = 1:tno_test    % spike number after first spike following 'BurstOn' event
    cspt = spt_test(k,st+1:st+nmbn);
    pki = find(cspt,1,'first');
    if isempty(pki)
        pspno_tt(k,1) = 0;
    else
        cspt = spt_test(k,pki+1:pki+nmbn);
        pspno_tt(k,1) = sum(cspt);
    end
end
spno = [sum(pspno) sum(pspno_tt)];
kn = st / nmbn + 1;
spno2 = spno2(1:kn);
spno2(kn) = spno(end);
if ~isequal(spno,spno2)
    disp('Technical error 122.')
    keyboard
end
% figure      % plot spike number
% plot(spno)

% Calculate p-value
p_value = makep(spno,kn);
keyboard

% -------------------------------------------------------------------------
function p_value = makep(spno,kn)
% Calculates p value for spike number.

nullhypkld = spno(1:kn-1);   % nullhypothesis
testkld = spno(kn);  % value to test
sno = length(nullhypkld(:));   % sample size for nullhyp. distribution
p_value = length(find(nullhypkld<=testkld)) / sno;