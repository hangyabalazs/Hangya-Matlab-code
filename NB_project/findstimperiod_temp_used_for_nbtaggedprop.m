function [activation_start activation_end activation_peak activation_time...
        baseline maxafter H] = findstimperiod(cellid,varargin)
%FINDSTIMPERIOD   Find time period of stimulated spikes.
%   [I1 I2] = FINDSTIMPERIOD(CELLID) seeks for an increase of firing rate
%   after an event. First, the program calculates Spike Densiti Function by
%   convolving the event-aligned raster plots with a variable Gaussian
%   window (see SOMPSTH_CALL2). Second, it finds maximal firing as maximum
%   SDF within a window from the event. Baseline firing is determined by
%   mean pre-event firing probability. Next, the time course of activation
%   is assessed by half-baseline crossings before and after the maximum.
%   This temporal window of activation is then returned in I1 (start) and
%   I2 (end timestamp).
%
%   [I1 I2 P] = FINDSTIMPERIOD(CELLID) returns the peak time of activation
%   (in seconds) as third output argument (P).
%
%   [I1 I2 P H] = FINDSTIMPERIOD(CELLID) plots the adaptive SDF and returns
%   the figure handle in H. Activation parameters are stored with the
%   figure as application data.
%
%   FINDSTIMPERIOD accepts the following parameter, value pairs as optional
%   input arguments (with default values):
%       'event', 'PulseOn' - the event to which the window is locked
%   	'window', [-0.005 0.01] - extent of baseline and test period
%           relative to the event in seconds; in seconds
%   	'margins', [-0.01 0.01] - margin in seconds used for extending the
%           'window' for adaptive SDF calculation to reduce convolution edge
%           effects; the SDF is later restricted to the 'window'
%       'dt', 0.0005 - time resolution of the bin raster; in seconds
%       'display', false - control of plotting event-locked raster plot
%       'valid_trials', 'all' - filter for trials; default: include all
%           trials; logical array or index set
%
%   See also FINDSEGS2 and ABS2RELTIMES.

%   Edit log: BH 4/25/12, 5/4/12

% Default arguments
default_args = {...
    'window',     [-0.005 0.01];... % window relative to the event, in seconds
    'margin',     [-0.01 0.01];...  % margins for psth calculation to get rid of edge effect due to smoothing
    'dt',         0.0005;...        % time resolution of the binraster, in seconds
    'display',    false;...         % control displaying rasters and PSTHs
    'event',      'PulseOn';...     % default event: 'PulseOn'
    'valid_trials','all';...        % valid trials - use all trials by default
    };
[g,error] = parse_args(default_args,varargin{:});

% Set parameters and load CellBase variables
ST = loadcb(cellid,'STIMSPIKES');   % load prealigned spikes for stimulation events
SE = loadcb(cellid,'StimEvents');
epoch_pos = findcellstr(ST.events(:,1),g.event);
if epoch_pos == 0
    error('Event not found.');
end
stimes = ST.event_stimes{epoch_pos};
time = (g.window(1)+g.margin(1)):g.dt:(g.window(2)+g.margin(2));

% Valid trials
valid_trials = parseValidTrials(SE,g.event,g.valid_trials);

% Calculate bin raster
spt = stimes2binraster(stimes(valid_trials),time,g.dt);

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
    
    % Plot raster plot and PSTH for 'PulseOn'
    set(gcf,'renderer','painters')   % temporarily change renderer because OpenGL locks the plot which result an error in legend layout handling
    viewcell2b(cellid,'TriggerName',g.event,'SortEvent',SEvent,'ShowEvents',ShEvent,'ShowEventsColors',{ShEvColors},...
        'FigureNum',FNum,'eventtype','stim','window',g.window,'dt',g.dt,'sigma',sigma,'PSTHstd',PSTHstd,'Partitions',parts,...
        'EventMarkerWidth',0,'PlotZeroLine','on')
    pause(0.05)   % if reset the renderer two early, the same error occurs
    set(gcf,'renderer','opengl')   % reset renderer
end

% PSTH
dtt = g.dt * 1000;   % resolution of bin raster in ms
wn = g.window * 1000;   % window boundaries in ms
margin = g.margin * 1000;   % margin in ms
if nargout < 7
    [activation_start activation_end activation_peak activation_time...
        baseline maxafter] = newpsth(spt,spt,dtt,wn,margin); %#ok<*ASGLU>
else
    [activation_start activation_end activation_peak activation_time...
        baseline maxafter H] = newpsth(spt,spt,dtt,wn,margin);
    tt = regexprep(cellid,'_',' ');     % add title to the figure
    title(tt)
end
activation_start = activation_start / 1000;   % convert back to seconds
activation_end = activation_end / 1000;
activation_peak = activation_peak / 1000;

% -------------------------------------------------------------------------
function [activation_start activation_end activation_peak activation_time...
    baseline_prob maxafter H] = newpsth(spt_baseline,spt_test,dt,win,margin)

% Control display (plot only if the plot handle is requested)
if nargout < 7
    dsply = false;
else
    dsply = true;
end

% Trial number and epoch length
[tno_baseline tl] = size(spt_baseline);
[tno_test tl] = size(spt_test);

% Pre-stimulus time window to consider for null hypothesis
stm = abs(win(1)+margin(1)) / dt;

% Merged spike train
sptb = spt_baseline(:,1:stm);
[x0 allspks_baseline] = find(sptb);
ts_baseline = sort(allspks_baseline)';

sptt = spt_test(:,stm+1:end);
[x0 allspks_test] = find(sptt);
ts_test = stm + sort(allspks_test)';
ts = [ts_baseline ts_test];

% Calculate adaptive SDF with variable Gaussian Kernel
prob = [sum(sptb)/tno_baseline sum(sptt)/tno_test]  / dt;  % prob. if 1 ms bins; dt is measured in ms!
spno = length(ts);
agvd = zeros(1,tl);
for t = 1:spno
    spi = ts(t);
    tspt = zeros(1,tl);
    tspt(spi) = 1;
    wbh = gausswin(9,prob(spi)*50);   % kernel
    wbh = wbh / sum(wbh);
    if spi <= stm
        agvd = agvd + filtfilt(wbh,1,tspt) / tno_baseline;   % convolution from both directions
    else
        agvd = agvd + filtfilt(wbh,1,tspt) / tno_test;
    end
end
ppsth_aconv = agvd / dt * 1000;   % SDF
psth_aconv = ppsth_aconv((stm+1+win(1)/dt):(stm+1+win(2)/dt));
stn = abs(win(1)) / dt;

% Plot SDF
time = win(1):dt:win(end);
if dsply
    H = figure;
    plot(time,psth_aconv,'k')
    xlim([time(1) time(end)])
end

% Activation time
baseline_prob = mean(prob(1:stm)) * 1000;  % spikes/sec (was spikes/bin before)
nst = length(psth_aconv);   % changing this allows further restricting the window
maxafter = max(psth_aconv(stn+1:nst));
if maxafter < baseline_prob     % putative activation, if firing goes above baseline
    activation_start = NaN;
    activation_end = NaN;
    activation_peak = NaN;
    activation_time = 0;   % if firing does not go above baseline
else
    maxinx = stn + find(psth_aconv(stn+1:nst)==maxafter,1,'first');   % maximal firing
    thr = baseline_prob + (maxafter - baseline_prob) / 2;
    pas = valuecrossing(time(stn+1:maxinx),psth_aconv(stn+1:maxinx),thr,'up');
    pas_inx = valuecrossing(stn+1:maxinx,psth_aconv(stn+1:maxinx),thr,'up');
    if isempty(pas)
        pas = time(stn+1);
        pas_inx = stn + 1;
    end
    pas_inx = round(pas_inx(end));
    activation_start = pas(end);   % last crossing of one and a half-baseline probability before maximum
    pae = valuecrossing(time(maxinx:nst),psth_aconv(maxinx:nst),thr,'down');
    pae_inx = valuecrossing(maxinx:nst,psth_aconv(maxinx:nst),thr,'down');
    if isempty(pae)
        pae = time(nst);
        pae_inx = nst;
    end
    pae_inx = round(pae_inx(1));
    activation_end = pae(1);   % first crossing of one and a half-baseline probability after maximum
    activation_time = activation_end - activation_start;
    activation_peak = time(maxinx) - time(stn+1);    % peak time of activation
end

% Plot
if dsply
    figure(H)
    hold on
    plot(time(pas_inx:pae_inx),psth_aconv(pas_inx:pae_inx),'Color',[0.8 0 0],'LineWidth',2)
    y_lim = ylim;
    ylim([0 y_lim(2)])
    xlabel('time [ms]')
    ylabel('firing rate')
    setappdata(H,'activation_start',activation_start)   % store the variables with the figure
    setappdata(H,'activation_end',activation_end)
    setappdata(H,'activation_peak',activation_peak)
    setappdata(H,'activation_time',activation_time)
    setappdata(H,'baseline',baseline_prob)
    setappdata(H,'maxafter',maxafter)
end