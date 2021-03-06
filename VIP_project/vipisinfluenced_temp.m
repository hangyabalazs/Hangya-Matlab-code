function [Wpi inhibition_start inhibition_end inhibition_time ...
    Wpa activation_start activation_end activation_time H] = ...
    vipisinfluenced(cellid)
%SOMISINHIBITED2   Test for event-locked inhibition.
%   P = SOMISINHIBITED2(CELLID) tests whether spike number after an event
%   (e.g. 'BurstOn') decreases. It uses a before-event baseline for
%   null-hypothesis estimataion of spike number. First, the program
%   calculates Spike Densiti Function by convolving the raster plots with a
%   variable Gaussian window (see SOMPSTH_CALL2). Second, it finds minimal
%   firing as minimum SDF within 100 ms from the event. Baseline firing is
%   determined by mean pre-event firing probability. Next, the time course
%   of inhibition is assessed by half-baseline crossings before and after
%   the minimum. This temporal window of inhibition is then used to bin the
%   baseline raster. Spike counts for baseline and spike counts in the
%   previously determined inhibition window are compared using Mann-Whitney
%   U-test (p-value is returned).
%
%   SOMISINHIBITED2 takes a different bin raster (typically aligned to
%   'BurstOn' events) for baseline distribution than the the one for
%   testing against baseline (typically aligned to 'PulseOn' events). Thus,
%   SDF reflects a hibrid spike raster merged at the event of interest.
%
%   [P S E T H] = SOMISINHIBITED2(CELLID) returns start (S) and end (E)
%   point of inhibition, total inhibition time (T) and figure handle for
%   SDF plot (H).
%
%   See also SOMISSTIM3 and SOMPSTH_CALL2.

% issues to solve: maximal time of inhibition before second pulse

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
% minfreq = min([TE.BurstNPulse]);
% maxpow = max([TE.PulsePower]);
% inx = ~isnan(getfield(TE,EventName2)) & TE.BurstNPulse==minfreq & TE.PulsePower==maxpow;
inx = ~isnan(getfield(TE,EventName2));
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

% PSTH
dtt = dt * 1000;   % resolution of bin raster in ms
wn = win * 1000;   % window boundaries in ms
pct = 0.5;         % percentage for threshold; 0.5 corresponde to half-baseline for inhibition and one and a half-baseline for activation
[inhibition_start inhibition_end inhibition_time Wpi ...
    activation_start activation_end activation_time Wpa H] ...
    = newpsth(spt1,spt2,dtt,wn,pct);

% -------------------------------------------------------------------------
function [inhibition_start inhibition_end inhibition_time Wpi ...
    activation_start activation_end activation_time Wpa H] ...
    = newpsth(spt_baseline,spt_test,dt,win,pct)

% Trial number and epoch length
[tno_baseline tl] = size(spt_baseline);
[tno_test tl] = size(spt_test);

% Pre-stimulus time window to consider for null hypothesis
st = abs(win(1)) / dt;   % in ms

% Merged spike train
sptb = spt_baseline(:,1:st);
[x0 allspks_baseline] = find(sptb);
ts_baseline = sort(allspks_baseline)';

sptt = spt_test(:,st+1:end);
[x0 allspks_test] = find(sptt);
ts_test = st + sort(allspks_test)';
ts = [ts_baseline ts_test];

% Calculate adaptive SDF with variable Gaussian Kernel
prob = [sum(sptb)/tno_baseline sum(sptt)/tno_test] / dt;  % prob. if 1 ms bins; dt is measured in ms!
spno = length(ts);
agvd1 = zeros(1,tl);
agvd2 = zeros(1,tl);
for t = 1:spno
    spi = ts(t);
    tspt = zeros(1,tl);
    tspt(spi) = 1;
    if prob(spi) > 1
        keyboard
    end
    wbh = gausswin(9,prob(spi)*50);   % kernel
    wbh = wbh / sum(wbh);
    if spi <= st
        agvd1 = agvd1 + filtfilt(wbh,1,tspt) / tno_baseline;   % convolution from both directions
    else
        agvd1 = agvd1 + filtfilt(wbh,1,tspt) / tno_test;
    end
end
psth_aconv = [agvd1(1:st) agvd1(st+1:end)] / dt * 1000;   % SDF

% Plot SDF
H = figure;
time = win(1):dt:win(end);
plot(time,psth_aconv,'k')
xlim([time(1) time(end)])

% Inhibition time
baseline_prob = mean(prob(1:st)) / dt * 1000;  % spikes/sec (was spikes/bin before)
nst = (abs(win(1)) + 100) / dt;   % index for 100 ms after stim event
minafter = min(psth_aconv(st+1:nst));
if minafter >= (baseline_prob * pct)  % inhibition, if firing goes below half-baseline
    inhibition_start = NaN;
    inhibition_end = NaN;
    inhibition_time = 0;   % if firing does not go below baseline
    Wpi = NaN;
else
    mininx = st + find(psth_aconv(st+1:nst)==minafter,1,'first');   % minimal firing
    pis = valuecrossing(time(st+1:mininx),psth_aconv(st+1:mininx),baseline_prob*pct,'down');
    pis_inx = valuecrossing(st+1:mininx,psth_aconv(st+1:mininx),baseline_prob*pct,'down');
    if isempty(pis)
        pis = time(st+1);
        pis_inx = st + 1;
    end
    pis_inx = round(pis_inx(end));
    inhibition_start = pis(end);   % last crossing of half-baseline probability before minimum
    pie = valuecrossing(time(mininx:nst),psth_aconv(mininx:nst),baseline_prob*pct,'up');
    pie_inx = valuecrossing(mininx:nst,psth_aconv(mininx:nst),baseline_prob*pct,'up');
    if isempty(pie)
        pie = time(nst);
        pie_inx = nst;
    end
    pie_inx = round(pie_inx(1));
    inhibition_end = pie(1);   % first crossing of half-baseline probability after minimum
    inhibition_time = inhibition_end - inhibition_start;
    
    % Nullhypothesis distribution
    wns = pie_inx - pis_inx + 1;
    spnon = sptb(:,st-wns*floor((st-1)/wns):st-1);
    spnon2 = reshape(spnon',1,tno_baseline*wns*floor((st-1)/wns))';
    spnon3 = reshape(spnon2,wns,tno_baseline*floor((st-1)/wns));
    spno_null = sum(spnon3);
    
    % Test distribution
    spno_test = sum(sptt(:,pis_inx-st:pie_inx-st),2);
    
    % Mann-Whitney test
    [Wpi,Whi] = b_ranksum2(spno_null,spno_test,'alpha',0.01);
    if Whi
        clri = [0 153 255] / 256;
    else
        clri = [102 255 255] / 256;
    end
end

% Activation time
maxafter = max(psth_aconv(st+1:nst));
if maxafter <= ((1 + pct) * baseline_prob)  % activation, if firing goes above one and a half-baseline
    activation_start = NaN;
    activation_end = NaN;
    activation_time = 0;   % if firing does not go above baseline
    Wpa = NaN;
else
    maxinx = st + find(psth_aconv(st+1:nst)==maxafter,1,'first');   % maximal firing
    pas = valuecrossing(time(st+1:maxinx),psth_aconv(st+1:maxinx),(1+pct)*baseline_prob,'up');
    pas_inx = valuecrossing(st+1:maxinx,psth_aconv(st+1:maxinx),(1+pct)*baseline_prob,'up');
    if isempty(pas)
        pas = time(st+1);
        pas_inx = st + 1;
    end
    pas_inx = round(pas_inx(end));
    activation_start = pas(end);   % last crossing of one and a half-baseline probability before maximum
    pae = valuecrossing(time(maxinx:nst),psth_aconv(maxinx:nst),(1+pct)*baseline_prob,'down');
    pae_inx = valuecrossing(maxinx:nst,psth_aconv(maxinx:nst),(1+pct)*baseline_prob,'down');
    if isempty(pae)
        pae = time(nst);
        pae_inx = nst;
    end
    pae_inx = round(pae_inx(1));
    activation_end = pae(1);   % first crossing of one and a half-baseline probability after maximum
    activation_time = activation_end - activation_start;
    
    % Nullhypothesis distribution
    wns = pae_inx - pas_inx + 1;
    spnon = sptb(:,st-wns*floor((st-1)/wns):st-1);
    spnon2 = reshape(spnon',1,tno_baseline*wns*floor((st-1)/wns))';
    spnon3 = reshape(spnon2,wns,tno_baseline*floor((st-1)/wns));
    spno_null = sum(spnon3);
    
    % Test distribution
    spno_test = sum(sptt(:,pas_inx-st:pae_inx-st),2);
    
    % Mann-Whitney test
    [Wpa,Wha] = b_ranksum2(spno_null,spno_test,'alpha',0.01);
    if Wha
        clra = 'red';
    else
        clra = [255 102 0] / 256;
    end
end

% Plot
if exist('clri','var')
    hold on
    plot(time(pis_inx:pie_inx),psth_aconv(pis_inx:pie_inx),'Color',clri,'LineWidth',2)
    x_lim = xlim;
    y_lim = ylim;
    text(x_lim(1)+(x_lim(2)-x_lim(1))*0.6,y_lim(1)+(y_lim(2)-y_lim(1))*0.8,['{\itMW test}, p = ',num2str(Wpi)],'Color',clri);
end
if exist('clra','var')
    hold on
    plot(time(pas_inx:pae_inx),psth_aconv(pas_inx:pae_inx),'Color',clra,'LineWidth',2)
    x_lim = xlim;
    y_lim = ylim;
    text(x_lim(1)+(x_lim(2)-x_lim(1))*0.6,y_lim(1)+(y_lim(2)-y_lim(1))*0.7,['{\itMW test}, p = ',num2str(Wpa)],'Color',clra);
end