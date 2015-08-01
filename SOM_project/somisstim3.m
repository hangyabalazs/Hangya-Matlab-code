function [p_value Idiff p_value2 Idiff2 p_value3 Idiff3 p_value4 Idiff4] = ...
    somisstim3(cellid)
%SOMISSTIM3   Assessment of optic tagging.
%   [P1 I1 P2 I2 P3 I3 P4 I4] = SOMISSTIM3(CELLID) calculates information
%   distances and corresponding p values for light tagging with four
%   methods for the cell given in CELLID (see CellBase documentation).
%   Briefly, a spike raster is calculated with 1 ms resolution. The raster
%   is broken down to 10 ms time bins and spike latency distribution for
%   first spikes is computed within each bin. Pairwise information
%   divergence measures are calculated for the before-light distributions
%   to form a null-hypothesis distribution for distances. The distances of
%   the first after-light distribution (using 'BurstOn' events - see
%   CellBase documentation) from all before-light distributions is
%   calculated and the median of these values is tested against the
%   null-hypothesis distribution. Output arguments:
%       P1 and I1 correspond to symmetric Kullback-Leibler divergence; 
%       P2 and I2 correspond to a modified symmetric Kullback-Leibler
%       divergence, where no normalisation is performed after restricting
%       the distributions to their common support
%       P3 and I3 correspond to Kullback-Leibler distances from uniform
%       distributions (with preserving the bin of no spikes in each
%       distribution); note, that no pairwise distances are calculated in
%       this case
%       P4 and I4 correspond to a modified version of Jensen-Shannon
%       divergence (see Endres and Schindelin, 2003)
%
%   SOMISSTIM3 takes a different bin raster (typically aligned to 'BurstOn'
%   events) for baseline distribution than the the one for testing against
%   baseline (typically aligned to 'PulseOn' events).
%
%   Reference:
%   Endres DM, Schindelin JE (2003) A new metric for probability
%   distributions. IEEE Transactions on Information Theory 49:1858-1860.
%
%   See also SOMISSTIM2_CALL, STIMES2BINRASTER, KLDISTSYM and JSDIV.

% Input argument check
if nargin < 2
    win = [-0.6 0.6];  % time window for bin raster (for Hyun: 1.2)
    dt = 0.001;   % resolution of bin raster in s (for Hyun: 0.0005)
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
dtt = dt * 1000;   % resolution of bin raster in ms
wn = win * 1000;   % window boundaries in ms
[p_value Idiff p_value2 Idiff2 p_value3 Idiff3 p_value4 Idiff4] = ...
    isikldist(spt1,spt2,dtt,wn,res);

% -------------------------------------------------------------------------
function [p_value Idiff p_value2 Idiff2 p_value3 Idiff3 p_value4 Idiff4] = ...
    isikldist(spt_baseline,spt_test,dt,win,res)

% Trial number and epoch length
[tno tl] = size(spt_baseline);

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
kld = nan(kn,kn);   % pairwise symmetric KL-distance
kld_mod = nan(kn,kn);   % pairwise modified symmetric KL-distance, with no normalization after restricting to the intersect of supports
kld_uniform = nan(1,kn);   % distance from 'uniform' distribution (with first bin preserved); not pairwise!
jsd = nan(kn,kn);  % pairwise modified JS-divergence (which is a metric!)
for k1 = 1:kn
    D1 = nhlsi(:,k1);
    for k2 = k1+1:kn
        D2 = nhlsi(:,k2);
        kld(k1,k2) = KLdistsym(D1,D2);
        kld_mod(k1,k2) = KLdistsym_nonorm(D1,D2);
        jsd(k1,k2) = sqrt(JSdiv(D1,D2)*2);
    end
    uniD = [D1(1); ones(length(D1)-1,1)*(1-D1(1))/(length(D1)-1)];
    kld_uniform(1,k1) = KLdistsym(D1,uniD);
end
% figure    % plot KL-distance
% imagesc(kld)

% Calculate p-value and information difference
[p_value Idiff] = makep(kld,kn);
[p_value2 Idiff2] = makep(kld_mod,kn);
[p_value3 Idiff3] = makep2(kld_uniform,kn);
[p_value4 Idiff4] = makep(jsd,kn);
% keyboard

% -------------------------------------------------------------------------
function D = KLdistsym_nonorm(P,Q)
%KLDISTSYM_NONORM   Modidied symmetric Kullbach-Leibler distance.
%   D = KLDISTSYM_NONORM(P,Q) calculates the symmetric Kullbach-Leibler
%   distance of the two input distributions (see KLDISTSYM), but performs
%   no normalization after restricting the distributions to their common
%   support.
%
%   See also KLDISTSYM.

% Input argument check
error(nargchk(2,2,nargin))
if abs(sum(P(:))-1) > 0.00001 || abs(sum(Q(:))-1) > 0.00001
    error('Input arguments must be probability distributions.')
end
if ~isequal(size(P),size(Q))
    error('Input distributions must be of the same size.')
end

% KL-distance
P2 = P(P.*Q>0);     % restrict to the common support
Q2 = Q(P.*Q>0);
D = sum(P2.*log(P2./Q2)) + sum(Q2.*log(Q2./P2));

% Alternative way of computation:
% HPQ = -sum(P2.*log(Q2));      % cross-entropy
% HP = -sum(P2.*log(P2));       % entropy
% D = HPQ - HP;

% -------------------------------------------------------------------------
function [p_value Idiff] = makep(kld,kn)
% Calculates p value from distance matrix.

pnhk = kld(1:kn-1,1:kn-1);
nullhypkld = pnhk(~isnan(pnhk));   % nullhypothesis
testkld = median(kld(1:kn-1,kn));  % value to test
sno = length(nullhypkld(:));   % sample size for nullhyp. distribution
p_value = length(find(nullhypkld>=testkld)) / sno;
Idiff = testkld - median(nullhypkld);

% -------------------------------------------------------------------------
function [p_value Idiff] = makep2(kld,kn)
% Calculates p value from distance matrix.

pnhk = kld(1:kn-1);
nullhypkld = pnhk(~isnan(pnhk));   % nullhypothesis
testkld = kld(kn);  % value to test
sno = length(nullhypkld(:));   % sample size for nullhyp. distribution
p_value = length(find(nullhypkld>=testkld)) / sno;
Idiff = testkld - median(nullhypkld);