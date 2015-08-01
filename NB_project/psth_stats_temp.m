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
%   	'testwin', [-0.25 0.1] - limits of baseline and test window for
%           statistical testing, time relative to 0 in seconds
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
%   See also VIPISINFLUENCED3B_A1, ULTIMATE_PSTH, APSTH2 and DAPSTH.

%   Edit log: BH 8/12/12

% Default arguments
prs = inputParser;
addRequired(prs,'spt',@isnumeric)   % bin raster
addRequired(prs,'psth',@isnumeric)   % PSTH
addRequired(prs,'dt',@isnumeric)   % time resolution of the binraster and PSTH, in seconds
addRequired(prs,'win',@(s)isnumeric(s)&isequal(length(s),2))  % time window relative to the event, in seconds
addParamValue(prs,'testwin',[-0.25 0.1],@(s)isnumeric(s)&isequal(length(s),2))  % time window relative to the event for stat. testing, in seconds
addParamValue(prs,'relative_threshold',0.5,@(s)isnumeric(s)&s>0&s<1)   % threshold used to assess interval limits
addParamValue(prs,'display',false,@(s)islogical(s)|ismember(s,[0 1]))   % control displaying rasters and PSTHs
parse(prs,spt,psth,dt,win,varargin{:})
g = prs.Results;

% Window for testing the potential effect
WN = g.testwin / dt;   % two-component vector, baseline and testwindow
WNb = abs(WN(1));   % baseline window
WNt = WN(2);   % test window

% Trial number
tno = size(spt,1);

% Time vector
time = win(1):dt:win(2);

% Baseline time window
st = abs(win(1)) / dt;   % in ms

% Spiking probability
sptb = spt(:,1:st);
sptt = spt(:,st+1:end);
probb = sum(sptb) / tno;   % spiking prob. (strictly) before time 0

% Inhibition time
baseline_prob = mean(probb(1:st)) / dt;  % spikes/sec (was spikes/bin before)
nst = st + WNt;   % index for end of test window
minafter = min(psth(st+1:nst));
if minafter > baseline_prob     % putative inhibition, if it goes below baseline
    inhibition_start = NaN;
    inhibition_end = NaN;
    inhibition_peak = NaN;
    inhibition_time = 0;   % if firing does not go below baseline
    Wpi = NaN;
else
    mininx = st + find(psth(st+1:nst)==minafter,1,'first');   % minimal firing
    thr = baseline_prob - (baseline_prob - minafter) * g.relative_threshold;  % threshold is determined in proportion of peak-baseline distance
    pis = valuecrossing(time(st+1:mininx),psth(st+1:mininx),thr,'down');
    pis_inx = valuecrossing(st+1:mininx,psth(st+1:mininx),thr,'down');
    if isempty(pis)
        pis = time(st+1);
        pis_inx = st + 1;
    end
    pis_inx = round(pis_inx(end));
    inhibition_start = pis(end);   % last crossing of half-baseline probability before minimum
    pie = valuecrossing(time(mininx:nst),psth(mininx:nst),thr,'up');
    pie_inx = valuecrossing(mininx:nst,psth(mininx:nst),thr,'up');
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
    wnnm = floor((st-10)/WNb);   % last 10 values omitted from baseline because of drop due to smoothing (no contrib. from after 0) 
    psp = nan(wns,wnnm*tno);
    for k = 1:wnnm
        inx = st-k*WNb-9:st-(k-1)*WNb-10;
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
maxafter = max(psth(st+1:nst));
if maxafter < baseline_prob     % putative activation, if firing goes above baseline
    activation_start = NaN;
    activation_end = NaN;
    activation_peak = NaN;
    activation_time = 0;   % if firing does not go above baseline
    Wpa = NaN;
else
    maxinx = st + find(psth(st+1:nst)==maxafter,1,'first');   % maximal firing
    thr = baseline_prob + (maxafter - baseline_prob) * g.relative_threshold;  % threshold is determined in proportion of peak-baseline distance
    pas = valuecrossing(time(st+1:maxinx),psth(st+1:maxinx),thr,'up');
    pas_inx = valuecrossing(st+1:maxinx,psth(st+1:maxinx),thr,'up');
    if isempty(pas)
        pas = time(st+1);
        pas_inx = st + 1;
    end
    pas_inx = round(pas_inx(end));
    activation_start = pas(end);   % last crossing of one and a half-baseline probability before maximum
    pae = valuecrossing(time(maxinx:nst),psth(maxinx:nst),thr,'down');
    pae_inx = valuecrossing(maxinx:nst,psth(maxinx:nst),thr,'down');
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
    wnnm = floor((st-10)/WNb);   % last 10 values omitted from baseline because of drop due to smoothing (no contrib. from after 0) 
    psp = nan(wns,wnnm*tno);
    for k = 1:wnnm
        inx = st-k*WNb-9:st-(k-1)*WNb-10;
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