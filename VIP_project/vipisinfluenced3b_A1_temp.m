function [baseline Wpi minvalue inhibition_start inhibition_end inhibition_peak inhibition_time ...
    Wpa maxvalue activation_start activation_end activation_peak activation_time H Hr] = ...
    vipisinfluenced3b(cellid)
%VIPISINFLUENCED3B   Test for event-locked inhibition and activation.
%   VIPISINFLUENCED3B(CELLID) tests whether spike number after an event
%   ('PulseOn') decreases. In VIPISINFLUENCED3B, there's no stimulus-free
%   baseline (designed for non-burst stimulation protocols).
%
%   First, the program calculates Spike Densiti Function by convolving the
%   raster plots with a variable Gaussian window (see SOMPSTH_CALL2).
%   Second, it finds minimal/maximal firing as minimum/maximum SDF within
%   200 ms from the event. Baseline firing is determined by mean pre-event
%   firing probability. Next, the time course of inhibition/activation is
%   assessed by crossings of the half-distence between the extreme and the
%   baseline before and after the minimum/maximum. This temporal window of
%   inhibition/activation is then used to bin the baseline raster. Spike
%   counts for baseline and spike counts in the previously determined
%   inhibition window are compared using Mann-Whitney U-test (p-value is
%   returned).
%
%   [B PI MI SI EI PTI TI PA MA SA EA PTA TA H] = VIPISINFLUENCED3B(CELLID)
%   returns minimal (MI) and maximal (MA) value within 100 ms from the
%   event; significance of inhibition (PI) and activation (PA); start (SI),
%   end (EI) and peak (PTI) time point of inhibition, total inhibition time
%   (TI), start point (SA), end (EA) and peak (PTA) time point of
%   activation, total activation time (TA) and figure handle for SDF plot
%   (H).
%
%   See also VIPISINFLUENCED, VIPISINFLUENCED2, SOMISINHIBITED2, SOMISSTIM3
%   and SOMPSTH_CALL2.

% Input argument check
if nargin < 2
    win = [-0.42 0.3];  % time window for bin raster
    dt = 0.005;   % resolution of bin raster in s
    dsply = 0;   % repress display
end

% Set parameters and load CellBase variables
EventName1 = 'PulseOn';
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

% Hibrid raster
rr = round(rand(1,size(spt2,1))*(size(spt1,1)-1)+0.5);
bspt = spt1(rr,time<0);
tspt = spt2(:,time>=0);
hspt = [bspt tspt];
Hr = rasterplot(hspt);

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
[baseline ...
    minvalue inhibition_start inhibition_end inhibition_peak inhibition_time Wpi ...
    maxvalue activation_start activation_end activation_peak activation_time Wpa H] ...
    = newpsth(spt1,spt2,dtt,wn);

% -------------------------------------------------------------------------
function [baseline_prob ...
    minafter inhibition_start inhibition_end inhibition_peak inhibition_time Wpi ...
    maxafter activation_start activation_end activation_peak activation_time Wpa H] ...
    = newpsth(spt_baseline,spt_test,dt,win)

% Window for testing the potential effect
WN = 200 / dt;   % 200 ms

% Trial number and epoch length
[tno_baseline tl] = size(spt_baseline);
[tno_test tl] = size(spt_test);

% Pre-stimulus time window to consider for null hypothesis
st = abs(win(1)) / dt;   % in ms

% Merged spike train
sptb = spt_baseline(:,1:st);
[x0 allspks_baseline] = find(sptb);
ts_baseline = sort(allspks_baseline)';

lag = 100 / dt;   % 100 ms lag to provide appropriate junction at 'st'
sptt = spt_test(:,st+1-lag:end);
[x0 allspks_test] = find(sptt);
ts_test = st - lag + sort(allspks_test)';

% Calculate adaptive SDF with variable Gaussian Kernel
probb = sum(sptb) / tno_baseline  / dt;  % prob. if 1 ms bins; dt is measured in ms!
probt = [zeros(1,st-lag) sum(sptt)/tno_test/dt];
spnob = length(ts_baseline);
spnot = length(ts_test);
agvd1 = zeros(1,tl);
agvd2 = zeros(1,tl);
for t = 1:spnob
    spi = ts_baseline(t);
    tspt = zeros(1,tl);
    tspt(spi) = 1;
    if probb(spi) > 1
        keyboard
    end
    wbh = gausswin(9,probb(spi)*50);   % kernel
    wbh = wbh / sum(wbh);
    agvd1 = agvd1 + filtfilt(wbh,1,tspt) / tno_baseline;   % convolution from both directions
end
for t = 1:spnot
    spi = ts_test(t);
    tspt = zeros(1,tl);
    tspt(spi) = 1;
    if probt(spi) > 1
        keyboard
    end
    wbh = gausswin(9,probt(spi)*50);   % kernel
    wbh = wbh / sum(wbh);
    agvd2 = agvd2 + filtfilt(wbh,1,tspt) / tno_test;   % convolution from both directions
end
psth_aconv = [agvd1(1:st) agvd2(st+1:end)] / dt * 1000;   % SDF

% Plot SDF
H = figure;
time = win(1):dt:win(end);
plot(time,psth_aconv,'k')
xlim([time(1) time(end)])

% Inhibition time
baseline_prob = mean(probb(1:st)) / dt * 1000;  % spikes/sec (was spikes/bin before)
nst = abs(win(1)) / dt + WN;   % index for 100 ms after stim event
minafter = min(psth_aconv(st+1:nst));
if minafter > baseline_prob     % putative inhibition, if it goes below baseline
    inhibition_start = NaN;
    inhibition_end = NaN;
    inhibition_peak = NaN;
    inhibition_time = 0;   % if firing does not go below baseline
    Wpi = NaN;
else
    mininx = st + find(psth_aconv(st+1:nst)==minafter,1,'first');   % minimal firing
    thr = baseline_prob - (baseline_prob - minafter) / 2;
    pis = valuecrossing(time(st+1:mininx),psth_aconv(st+1:mininx),thr,'down');
    pis_inx = valuecrossing(st+1:mininx,psth_aconv(st+1:mininx),thr,'down');
    if isempty(pis)
        pis = time(st+1);
        pis_inx = st + 1;
    end
    pis_inx = round(pis_inx(end));
    inhibition_start = pis(end);   % last crossing of half-baseline probability before minimum
    pie = valuecrossing(time(mininx:nst),psth_aconv(mininx:nst),thr,'up');
    pie_inx = valuecrossing(mininx:nst,psth_aconv(mininx:nst),thr,'up');
    if isempty(pie)
        pie = time(nst);
        pie_inx = nst;
    end
    pie_inx = round(pie_inx(1));
    inhibition_end = pie(1);   % first crossing of half-baseline probability after minimum
    inhibition_time = inhibition_end - inhibition_start;
    inhibition_peak = time(mininx) - time(st+1);    % peak time of inhibition
    
    % Nullhypothesis distribution
    wns = pie_inx - pis_inx + 1;
    wnnm = floor((st-10)/WN);   % last 10 values omitted from baseline because of drop due to smoothing (no contrib. from after 0) 
    psp = nan(wns,wnnm*tno_baseline);
    for k = 1:wnnm
        inx = st-k*WN-9:st-(k-1)*WN-10;
        cwn = sptb(:,inx);
        cwnps = psth_aconv(inx);
        mcw = find(cwnps==min(cwnps));
        mcw = mcw(1);
        inx2 = mcw-floor(wns/2):mcw+ceil(wns/2)-1;
        inx2 = inx2 - min(0,inx2(1)-1) - (max(length(inx),inx2(end)) - length(inx));
        if ~ismember(mcw,inx2)
            error('vipisinfluenced2:nullhypoIndexing','Programming error.')
        end
        mcwn = cwn(:,inx2);
        if ~isequal(size(mcwn,2),wns)
            error('vipisinfluenced2:nullhypoIndexing','Programming error.')
        end
        psp(:,(k-1)*tno_baseline+1:k*tno_baseline) = mcwn';
    end
    if any(isnan(psp))
        error('vipisinfluenced2:nullhypoIndexing','Programming error.')
    end
    spno_null = sum(psp);
    
    % Test distribution
    spno_test = sum(sptt(:,pis_inx-st+lag:pie_inx-st+lag),2);
    
    % Mann-Whitney test
    [Wpi,Whi] = b_ranksum2(spno_null,spno_test,'alpha',0.01);
    ranks = tiedrank([spno_null spno_test']);
    tp = length(spno_null);
    nullranksum = mean(ranks(1:tp));
    testranksum = mean(ranks(tp+1:end));
    if testranksum > nullranksum    % one-sided test
        Wpi = NaN;
        Whi = 0;
    end
    if Whi
        clri = [0 153 255] / 256;
    else
        clri = [102 255 255] / 256;
    end
end

% Activation time
maxafter = max(psth_aconv(st+1:nst));
if maxafter < baseline_prob     % putative activation, if firing goes above baseline
    activation_start = NaN;
    activation_end = NaN;
    activation_peak = NaN;
    activation_time = 0;   % if firing does not go above baseline
    Wpa = NaN;
else
    maxinx = st + find(psth_aconv(st+1:nst)==maxafter,1,'first');   % maximal firing
    thr = baseline_prob + (maxafter - baseline_prob) / 2;
    pas = valuecrossing(time(st+1:maxinx),psth_aconv(st+1:maxinx),thr,'up');
    pas_inx = valuecrossing(st+1:maxinx,psth_aconv(st+1:maxinx),thr,'up');
    if isempty(pas)
        pas = time(st+1);
        pas_inx = st + 1;
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
    activation_peak = time(maxinx) - time(st+1);    % peak time of activation
    
    % Nullhypothesis distribution
    wns = pae_inx - pas_inx + 1;
    wnnm = floor((st-10)/WN);   % last 10 values omitted from baseline because of drop due to smoothing (no contrib. from after 0) 
    psp = nan(wns,wnnm*tno_baseline);
    for k = 1:wnnm
        inx = st-k*WN-9:st-(k-1)*WN-10;
        cwn = sptb(:,inx);
        cwnps = psth_aconv(inx);
        mcw = find(cwnps==max(cwnps));
        mcw = mcw(1);
        inx2 = mcw-floor(wns/2):mcw+ceil(wns/2)-1;
        inx2 = inx2 - min(0,inx2(1)-1) - (max(length(inx),inx2(end)) - length(inx));
        if ~ismember(mcw,inx2)
            error('vipisinfluenced2:nullhypoIndexing','Programming error.')
        end
        mcwn = cwn(:,inx2);
        if ~isequal(size(mcwn,2),wns)
            error('vipisinfluenced2:nullhypoIndexing','Programming error.')
        end
        psp(:,(k-1)*tno_baseline+1:k*tno_baseline) = mcwn';
    end
    if any(isnan(psp))
        error('vipisinfluenced2:nullhypoIndexing','Programming error.')
    end
    spno_null = sum(psp);
    
    % Test distribution
    spno_test = sum(sptt(:,pas_inx-st+lag:pae_inx-st+lag),2);
    
    % Mann-Whitney test
    [Wpa,Wha] = b_ranksum2(spno_null,spno_test,'alpha',0.01);
    ranks = tiedrank([spno_null spno_test']);
    tp = length(spno_null);
    nullranksum = mean(ranks(1:tp));
    testranksum = mean(ranks(tp+1:end));
    if testranksum < nullranksum    % one-sided test
        Wpa = NaN;
        Wha = 0;
    end
    if Wha
        clra = 'red';
    else
        clra = [255 102 0] / 256;
    end
end

% Plot
figure(H)
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