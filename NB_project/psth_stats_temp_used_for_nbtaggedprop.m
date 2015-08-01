function stats = psth_stats(spt,psth,dt,win,varargin)
%PSTH_STATS   Statistics for firing rate change.
%   PSTH_STATS tests whether spike number after an event changes.
%
%   STATS = PSTH_STATS(SPT,PSTH,DT,WIN) finds minimal/maximal firing as
%   minimum/maximum PSTH within 100 ms from time zero (determined by time
%   limits of the window around 0, WIN and temporal resolution, DT).
%   Baseline firing is determined by mean firing probability from -250 ms 
%   to 0 using SPT, the bin raster corrsponding to the PSTH. Next, the time
%   course of inhibition/activation is assessed by crossings of the
%   half-distence between the extreme and the baseline before and after the
%   minimum/maximum. This temporal window of inhibition/activation is then
%   used to find corresponding intervals around local extremes in the
%   baseline raster. Spike counts for baseline and spike counts in the
%   previously determined inhibition window are compared using Mann-Whitney
%   U-test (p-value is calculated).
%
%   Default behavior of PSTH_STATS can be modified by using a set of
%   paramter-value pairs as optional input parameters. The following
%   parameters are implemented (with default values):
%   	'baselinewin', [-0.25 0] - limits of baseline window for 
%           statistical testing, time relative to 0 in seconds
%   	'testwin', [0 0.1] - limits of test window for statistical testing,
%           time relative to 0 in seconds
%       'relative_threshold', 0.5 - threshold used to assess start and end
%           points of activation and inhibition intervals; in proportion of
%           the peak-baseline difference
%       'display', false - controls plotting.
%
%   The output structure STATS contains the following fields:
%       'baseline' - baseline firing probability
%       'minvalue' - minimal firing rate in the test window
%       'inhibition_start' - start time of inhibition in seconds
%       'inhibition_end' - end time of inhibition in seconds
%       'inhibition_peak' - peak time of inhibition in seconds
%       'inhibition_time' - duration of inhibition in seconds
%       'Wpi' - p value for Mann-Whitney test for significant inhibition
%       'maxvalue' - maximal firing rate in the test window
%       'activation_start' - start time of activation in seconds
%       'activation_end' - end time of activation in seconds
%       'activation_peak' - peak time of activation in seconds
%       'activation_time' - duration of activation in seconds
%       'Wpa' - p value for Mann-Whitney test for significant activation
%
%   IMPORTANT NOTE: Please note that the baseline window is split to
%   smaller windows of the size of the test window. Use baseline windows
%   that are somewhat bigger than an integer multiple of the test window
%   length. The baseline window will be cropped BEFORE the nearest integer
%   multiple of the test window size. For instance, a test window of [0
%   0.1] and a baseline window of [-0.22 0] will result in an effective
%   baseline window of [-0.2 0], corresponding to twice the size of the
%   baseline window.
%
%   See also VIPISINFLUENCED3B_A1, ULTIMATE_PSTH, APSTH2 and DAPSTH.

%   Edit log: BH 8/12/12, 8/27/12, 7/15/13

% Default arguments
prs = inputParser;
addRequired(prs,'spt',@isnumeric)   % bin raster
addRequired(prs,'psth',@isnumeric)   % PSTH
addRequired(prs,'dt',@isnumeric)   % time resolution of the binraster and PSTH, in seconds
addRequired(prs,'win',@(s)isnumeric(s)&isequal(length(s),2))  % time window relative to the event, in seconds
addParamValue(prs,'baselinewin',[-0.25 0],@(s)isnumeric(s)&isequal(length(s),2))  % baseline time window relative to the event for stat. testing, in seconds
addParamValue(prs,'testwin',[0 0.1],@(s)isnumeric(s)&isequal(length(s),2))  % test time window relative to the event for stat. testing, in seconds
addParamValue(prs,'relative_threshold',0.5,@(s)isnumeric(s)&s>0&s<1)   % threshold used to assess interval limits
addParamValue(prs,'display',false,@(s)islogical(s)|ismember(s,[0 1]))   % control displaying rasters and PSTHs
parse(prs,spt,psth,dt,win,varargin{:})
g = prs.Results;

% Error handling
twlen = diff(g.testwin);   % length of the test window
bwlen = diff(g.baselinewin);   % length of the baseline window
if bwlen < twlen
    warning('psth_stats:inputArg','PSTH_STATS: For a valid MW-test, the baseline window should be longer than the test window.')
end
if bwlen > -g.win(1)
    error('psth_stats:inputArg','PSTH_STATS: The PSTH window does not cover the baseline window.')
end
if twlen > g.win(2)
    error('psth_stats:inputArg','PSTH_STATS: The PSTH window does not cover the test window.')
end

% Index for time 0
st = abs(win(1)) / dt;   % in ms
nullindex = st + 1;

% Window for testing the potential effect
WNb = g.baselinewin / dt + nullindex - 1;   % baseline window; convert to indices
WNt = g.testwin / dt + nullindex;   % test window; convert to indices
WNb = round(WNb);
WNt = round(WNt);
lWNb = WNb(2) - WNb(1) + 1;   % length of baseline window
lWNt = WNt(2) - WNt(1) + 1;   % length of test window

% Trial number
tno = size(spt,1);

% Time vector
time = win(1):dt:win(2);

% Spiking probability
sptb = spt(:,1:st);    % st < time 0
sptt = spt(:,st+1:end);   % st+1 = time 0
probb = sum(sptb) / tno;   % spiking prob. (strictly) before time 0

% Inhibition time
baseline_prob = mean(probb(WNb(1):WNb(2))) / dt;  % spikes/sec (was spikes/bin before)
minafter = min(psth(WNt(1):WNt(2)));
if minafter > baseline_prob     % putative inhibition, if it goes below baseline
    inhibition_start = NaN;
    inhibition_end = NaN;
    inhibition_peak = NaN;
    inhibition_time = 0;   % if firing does not go below baseline
    Wpi = NaN;
else
    mininx = st + find(psth(WNt(1):WNt(2))==minafter,1,'first');   % minimal firing
    thr = baseline_prob - (baseline_prob - minafter) * g.relative_threshold;  % threshold is determined in proportion of peak-baseline distance
    pis = valuecrossing(time(WNt(1):mininx),psth(WNt(1):mininx),thr,'down');
    pis_inx = valuecrossing(WNt(1):mininx,psth(WNt(1):mininx),thr,'down');
    if isempty(pis)
        pis = time(WNt(1));
        pis_inx = WNt(1);
    end
    pis_inx = round(pis_inx(end));
    inhibition_start = pis(end);   % last crossing of half-baseline probability before minimum
    pie = valuecrossing(time(mininx:WNt(2)),psth(mininx:WNt(2)),thr,'up');
    pie_inx = valuecrossing(mininx:WNt(2),psth(mininx:WNt(2)),thr,'up');
    if isempty(pie)
        pie = time(WNt(2));
        pie_inx = WNt(2);
    end
    pie_inx = round(pie_inx(1));
    inhibition_end = pie(1);   % first crossing of half-baseline probability after minimum
    inhibition_time = inhibition_end - inhibition_start;
    inhibition_peak = time(mininx) - time(st+1);    % peak time of inhibition
    
    % Nullhypothesis distribution
    wns = pie_inx - pis_inx + 1;
    wnnm = floor(lWNt/lWNb);   % split up the baseline according to test window length
    psp = nan(wns,wnnm*tno);
    for k = 1:wnnm
        inx = WNb(2)-k*lWNt+1:WNb(2)-(k-1)*lWNt;
        cwn = sptb(:,inx);
        cwnps = psth(inx);
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
        psp(:,(k-1)*tno+1:k*tno) = mcwn';
    end
    if any(isnan(psp))
        error('vipisinfluenced2:nullhypoIndexing','Programming error.')
    end
    spno_null = sum(psp);
    
    % Test distribution
    spno_test = sum(sptt(:,pis_inx-st:pie_inx-st),2);
    
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
maxafter = max(psth(WNt(1):WNt(2)));
if maxafter < baseline_prob     % putative activation, if firing goes above baseline
    activation_start = NaN;
    activation_end = NaN;
    activation_peak = NaN;
    activation_time = 0;   % if firing does not go above baseline
    Wpa = NaN;
else
    maxinx = st + find(psth(WNt(1):WNt(2))==maxafter,1,'first');   % maximal firing
    thr = baseline_prob + (maxafter - baseline_prob) * g.relative_threshold;  % threshold is determined in proportion of peak-baseline distance
    pas = valuecrossing(time(WNt(1):maxinx),psth(WNt(1):maxinx),thr,'up');
    pas_inx = valuecrossing(WNt(1):maxinx,psth(WNt(1):maxinx),thr,'up');
    if isempty(pas)
        pas = time(WNt(1));
        pas_inx = WNt(1);
    end
    pas_inx = round(pas_inx(end));
    activation_start = pas(end);   % last crossing of one and a half-baseline probability before maximum
    pae = valuecrossing(time(maxinx:WNt(2)),psth(maxinx:WNt(2)),thr,'down');
    pae_inx = valuecrossing(maxinx:WNt(2),psth(maxinx:WNt(2)),thr,'down');
    if isempty(pae)
        pae = time(WNt(2));
        pae_inx = WNt(2);
    end
    pae_inx = round(pae_inx(1));
    activation_end = pae(1);   % first crossing of one and a half-baseline probability after maximum
    activation_time = activation_end - activation_start;
    activation_peak = time(maxinx) - time(st+1);    % peak time of activation
    
    % Nullhypothesis distribution
    wns = pae_inx - pas_inx + 1;
    wnnm = floor(lWNb/lWNt);   % split up the baseline according to test window length
    psp = nan(wns,wnnm*tno);
    for k = 1:wnnm
        inx = WNb(2)-k*lWNt+1:WNb(2)-(k-1)*lWNt;
        cwn = sptb(:,inx);
        cwnps = psth(inx);
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
        psp(:,(k-1)*tno+1:k*tno) = mcwn';
    end
    if any(isnan(psp))
        error('vipisinfluenced2:nullhypoIndexing','Programming error.')
    end
    spno_null = sum(psp);
    
    % Test distribution
    spno_test = sum(sptt(:,pas_inx-st:pae_inx-st),2);
    
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
if g.display
    figure;
    plot(time,psth,'k')
    xlim([time(1) time(end)])
    if exist('clri','var')
        hold on
        plot(time(pis_inx:pie_inx),psth(pis_inx:pie_inx),'Color',clri,'LineWidth',2)
        x_lim = xlim;
        y_lim = ylim;
        text(x_lim(1)+(x_lim(2)-x_lim(1))*0.6,y_lim(1)+(y_lim(2)-y_lim(1))*0.8,['{\itMW test}, p = ',num2str(Wpi)],'Color',clri);
    end
    if exist('clra','var')
        hold on
        plot(time(pas_inx:pae_inx),psth(pas_inx:pae_inx),'Color',clra,'LineWidth',2)
        x_lim = xlim;
        y_lim = ylim;
        text(x_lim(1)+(x_lim(2)-x_lim(1))*0.6,y_lim(1)+(y_lim(2)-y_lim(1))*0.7,['{\itMW test}, p = ',num2str(Wpa)],'Color',clra);
    end
end

% Output statistics
stats.baseline = baseline_prob;
stats.minvalue = minafter;
stats.inhibition_start = inhibition_start;
stats.inhibition_end = inhibition_end;
stats.inhibition_peak = inhibition_peak;
stats.inhibition_time = inhibition_time;
stats.Wpi = Wpi;
stats.maxvalue = maxafter;
stats.activation_start = activation_start;
stats.activation_end = activation_end;
stats.activation_peak = activation_peak;
stats.activation_time = activation_time;
stats.Wpa = Wpa;