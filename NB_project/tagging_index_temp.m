function [p_value Idiff] = tagging_index(cellid,varargin)
%TAGGING_INDEX   Assessment of optic tagging.
%   [P I] = TAGGING_INDEX(CELLID) calculates information distances and
%   corresponding p values for light tagging for the cell given in CELLID
%   (see CellBase documentation). Briefly, a spike raster is calculated
%   with 1 ms resolution. The raster is devided to 10 ms time bins and
%   spike latency distribution for first spikes is computed within each
%   bin. Pairwise information divergence measures are calculated for the
%   before-light distributions to form a null-hypothesis distribution for
%   distances. The distances of the first after-light distribution (using
%   'BurstOn' events - see CellBase documentation) from all before-light
%   distributions is calculated and the median of these values is tested
%   against the null-hypothesis distribution. Output arguments P and I
%   correspond to a modified version of Jensen-Shannon divergence (see
%   Endres and Schindelin, 2003).
%
%   TAGGIN_INDEX takes a different bin raster (typically aligned to
%   'BurstOn' events) for baseline distribution than the the one for
%   testing against baseline (typically aligned to 'PulseOn' events).
%   Number of test trials used is maximized to 5000 by default.
%
%   Default behavior of TAGGING_INDEX can be modified by using a set of
%   paramter, value pairs as optional input parameters. The folowing
%   parameters are implemented (with default values):
%       'event', 'PulseOn' - the event to which the window is locked
%   	'window', [-0.6 0.6] - extent of baseline and test period
%           relative to the event in seconds; in seconds
%   	'dt', 0.001 - time resolution of the bin raster; in seconds
%       'display', false - control of plotting event-locked raster plot
%       'event_filter', 'none' - filter light-stimulation trials; see
%           implemented filter types below
%       'maxtrialno', 5000 - maximal number of light-stimulation trials
%           included; if ther are more valid trials, they are randomly
%           down-sampled
%
%   Implemented event filters:
%       'BurstNPulse' - only light pulse trains with a specified number of
%           pulses are used; this number has to be specified as a
%           'BurstNPulse', N (numeric) parameter, value pair
%       'BurstNPulse_maxPower' - the same as 'BurstNPulse' but light bursts
%           are further restricted to pulses with maximal power applied
%       'minNPulse_maxPower' - the light bursts with minimal number of
%           pulses and maximal power are used
%
%   Reference:
%   Endres DM, Schindelin JE (2003) A new metric for probability
%   distributions. IEEE Transactions on Information Theory 49:1858-1860.
%
%   See also NBTAGGING, SOMISSTIM3, STIMES2BINRASTER and JSDIV.


% Default arguments
default_args = {...
    'window',     [-0.6 0.6];... % time window for bin raster relative to the event, in seconds
    'dt',         0.001;...      % time resolution of the binraster, in seconds
    'display',    false;...      % control displaying rasters and PSTHs
    'event',      'PulseOn';...  % default reference event: 'PulseOn'
    'event_filter','none';...    % filter events based on properties
    'maxtrialno', 5000;...       % downsample events if more than 'maxtrialno'
    };
[g,error] = parse_args(default_args,varargin{:});

% Set parameters and load CellBase variables
EventName1 = 'BurstOn';     % for baseline
EventName2 = g.event;
ST = loadcb(cellid,'STIMSPIKES');   % load prealigned spikes for stimulation events
TE = loadcb(cellid,'StimEvents');

% Spike times for baseline period
epoch_pos1 = findcellstr(ST.events(:,1),EventName1);
epoch_pos2 = findcellstr(ST.events(:,1),EventName2);
if epoch_pos1 == 0 || epoch_pos2 == 0
    error('Epoch name not found');
end
stimes1 = ST.event_stimes{epoch_pos1};
stimes2 = ST.event_stimes{epoch_pos2};
time = g.window(1):g.dt:g.window(end);
valid_trials1 = find(~isnan(TE.(EventName1)));

% Spike times for test period - filter events
switch g.event_filter
    case 'minNPulse_maxPower'
        minfreq = min([TE.BurstNPulse]);
        maxpow = max([TE.PulsePower]);
        inx = ~isnan(TE.(EventName2)) & TE.BurstNPulse==minfreq & TE.PulsePower==maxpow;
        valid_trials2 = find(inx);
    case 'BurstNPulse_maxPower'
        maxpow = max([TE.PulsePower]);
        inx = ~isnan(TE.(EventName2)) & TE.BurstNPulse==g.BurstNPulse & TE.PulsePower==maxpow;
        valid_trials2 = find(inx);
    case 'BurstNPulse'
        inx = ~isnan(TE.(EventName2)) & TE.BurstNPulse==g.BurstNPulse;
        valid_trials2 = find(inx);
    case 'none'
        valid_trials2 = find(~isnan(TE.(EventName2)));
    otherwise
        error('Unknown event filter.')
end

% Downsaple if too many pulses
lm = g.maxtrialno;
if length(valid_trials2) > lm
    rp = randperm(length(valid_trials2));
    valid_trials2 = valid_trials2(sort(rp(1:lm)));
end

% Calculate bin rasters
spt1 = stimes2binraster(stimes1(valid_trials1),time,g.dt); %#ok<FNDSB>
spt2 = stimes2binraster(stimes2(valid_trials2),time,g.dt);

% Set input arguments for rater plot and PSTH
if g.display
    SEvent = 'BurstOff';
    FNum = 2;
    parts = 'all';
    sigma = 0.001;
    PSTHstd = 'on';
    ShEvent = {{'BurstOff'}};
    ShEvColors = hsv(length(ShEvent{1}));
    ShEvColors = mat2cell(ShEvColors,ones(size(ShEvColors,1),1),3);
    
    % Plot raster plot and PSTH for 'BurstOn'
    figure
    set(gcf,'renderer','painters')   % temporaray change renderer because OpenGL locks the plot which result an error in legend layout handling
    viewcell2b(cellid,'TriggerName',EventName1,'SortEvent',SEvent,'ShowEvents',ShEvent,'ShowEventsColors',{ShEvColors},...
        'FigureNum',FNum,'eventtype','stim','window',win,'dt',dt,'sigma',sigma,'PSTHstd',PSTHstd,'Partitions',parts,...
        'EventMarkerWidth',0,'PlotZeroLine','on')
    pause(0.05)   % if reset the renderer two early, the same error occurs
    set(gcf,'renderer','opengl')   % reset renderer
    
    % Plot raster plot and PSTH for 'PulseOn'
    figure
    set(gcf,'renderer','painters')   % temporaray change renderer because OpenGL locks the plot which result an error in legend layout handling
    viewcell2b(cellid,'TriggerName',EventName2,'SortEvent',SEvent,'ShowEvents',ShEvent,'ShowEventsColors',{ShEvColors},...
        'FigureNum',FNum,'eventtype','stim','window',win,'dt',dt,'sigma',sigma,'PSTHstd',PSTHstd,'Partitions',parts,...
        'EventMarkerWidth',0,'PlotZeroLine','on')
    pause(0.05)   % if reset the renderer two early, the same error occurs
    set(gcf,'renderer','opengl')   % reset renderer
end

% Calculate information distances and p values
res = 10;   % resolution in ms
dtt = g.dt * 1000;   % resolution of bin raster in ms
wn = g.window * 1000;   % window boundaries in ms
[p_value Idiff] = isikldist(spt1,spt2,dtt,wn,res);

% -------------------------------------------------------------------------
function [p_value Idiff] = isikldist(spt_baseline,spt_test,dt,win,res)

% Trial number and epoch length
[tno tl] = size(spt_baseline); %#ok<NASGU>

% Number of bins for ISI histograms
nmbn = round(res/dt);

% Pre-stimulus time window to consider for null hypothesis
st = abs(win(1)) / dt;   % number of pre-stim values in 'spt'

% ISI histogram - baseline
edges = 0:nmbn+1;
nm = floor(st/nmbn);
lsi = zeros(tno,nm);   % ISI's
slsi = zeros(tno,nm);  % sorted ISI's
hlsi = zeros(nmbn+1,nm);    % ISI hist.; +1: zero when no spike in the segment
nhlsi = zeros(nmbn+1,nm);   % normalized ISI histogram 
next = 1;
for t = 1:nmbn:st
    for k = 1:tno
        cspt = spt_baseline(k,t:t+nmbn-1);
        pki = find(cspt,1,'first');
        if ~isempty(pki)
            lsi(k,next) = pki;
        else
            lsi(k,next) = 0;
        end
    end
    slsi(:,next) = sort(lsi(:,next));
    hst = hist(slsi(:,next),edges);
    hlsi(:,next) = hst(1:end-1);
    nhlsi(:,next) = hlsi(:,next) / sum(hlsi(:,next));
    next = next + 1;
end

% ISI histogram - test
tno_test = size(spt_test,1);
lsi_tt = nan(tno_test,1);
for k = 1:tno_test
    cspt = spt_test(k,st+1:st+nmbn);
    pki = find(cspt,1,'first');
    if ~isempty(pki)
        lsi_tt(k,1) = pki;
    else
        lsi_tt(k,1) = 0;
    end
end
slsi_tt = sort(lsi_tt(:,1));
hst = hist(slsi_tt,edges);
hlsi(:,next) = hst(1:end-1);
nhlsi(:,next) = hlsi(:,next) / sum(hlsi(:,next));

% figure      % plot ISIs
% imagesc(lsi)
% figure      % plot sorted ISIs
% imagesc(slsi)
% figure      % plot ISI histograms
% imagesc(hlsi(2:end,:))

% Symmetric KL-divergence and JS-divergence
kn = st / nmbn + 1;
jsd = nan(kn,kn);  % pairwise modified JS-divergence (which is a metric!)
for k1 = 1:kn
    D1 = nhlsi(:,k1);
    for k2 = k1+1:kn
        D2 = nhlsi(:,k2);
        jsd(k1,k2) = sqrt(JSdiv(D1,D2)*2);
    end
end
% figure    % plot KL-distance
% imagesc(kld)

% Calculate p-value and information difference
[p_value Idiff] = makep(jsd,kn);
% keyboard

% -------------------------------------------------------------------------
function [p_value Idiff] = makep(kld,kn)
% Calculates p value from distance matrix.

pnhk = kld(1:kn-1,1:kn-1);
nullhypkld = pnhk(~isnan(pnhk));   % nullhypothesis
testkld = median(kld(1:kn-1,kn));  % value to test
sno = length(nullhypkld(:));   % sample size for nullhyp. distribution
p_value = length(find(nullhypkld>=testkld)) / sno;
Idiff = testkld - median(nullhypkld);