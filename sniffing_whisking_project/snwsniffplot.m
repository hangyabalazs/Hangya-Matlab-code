function snwsniffplot(rrat,rdate,csc1,csc2,msgn,thr,wl)
%SNWSNIFFPLOT   Whisking phase relative to respiration.
%   SNWSNIFFPLOT is a modified version of SNWPHASE4. It generates a
%   strucutre calles 'segments', which saves information about each segment
%   of a particular type of phase locking (see SNWPHASE4) for later
%   plotting (see SNWSEGMENTPLOTTER).
%
%   SNWSNIFFPLOT(RAT,DATE,CSC1,CSC2,MSGN,THR,WL) accepts input arguments
%   determining the animal (RAT) and session (DATE) information, which data
%   channels are used (CSC1 and CSC2 for sniffing and whisking,
%   respectively), a scalar determining whether sniffing is inverted
%   (MSGN=-1) or not (MSGN=1), threshold for RMS discrimination (THR) and
%   window size for RMS in ms (WL).
%
%   See also SNWPHASE4 and SNWSEGMENTPLOTTER.

% Input arguments
error(nargchk(0,7,nargin))
if nargin < 1
    rrat = 'R7';    % animal
end
if nargin < 2
    rdate = '2008-09-03_15-59-08';   % session
end
if nargin < 3
    csc1 = 'CSC4.ncs';  % respiration (sniffing) channel
end
if nargin < 4
    csc2 = 'CSC6.ncs';  % EMG (whisking) channel
end
if nargin < 5
    msgn = -1;   % control the polarity of the sniffing signal
end
if nargin < 6
    thr = 550;  % threshold for RMS discrimination
end
if nargin < 5
    wl = 20;    % 10 ms windows for RMS
end

% Define results directory
global DATAPATH
global DATADIR
resdir = [DATAPATH 'SniffWhisk\phasehist8\'];

% Load discriminated raw whisking
fnme = [DATAPATH 'SniffWhisk\disc\DISC_' rrat '_' rdate '.mat'];
load(fnme)
vdisc_rawwhisk = vdisc_rawwhisk / wl;

% Data import
plotting = false;
dsr = 1000;   % downsample rate

fn = ([DATADIR 'HSW\' rrat '\' rdate '\']);
[TimeStamp, ChanNum, SampleFrequency, NumValSamples, Samples, NlxHeader] = ...
    Nlx2MatCSC([fn csc1],[1 1 1 1 1],1,1,1);  % sniffing
sniff = msgn * Samples(:);
sr_sniff = 512 / mean(diff(TimeStamp)) * 10^6;    % original sampling rate (changed later)
dt_sniff = 1 / sr_sniff;
starttime1 = TimeStamp(1) / 10^6;
time_sniff = (0:length(sniff)-1) * dt_sniff + starttime1;
[TimeStamp, ChanNum, SampleFrequency, NumValSamples, Samples, NlxHeader] = ...
    Nlx2MatCSC([fn csc2],[1 1 1 1 1],1,1,1);  % whisking
whisk = Samples(:);
sr_whisk = 512 / mean(diff(TimeStamp)) * 10^6;    % original sampling rate for EMG (different from that for respiration in R7)
dt_whisk = 1 / sr_whisk;
starttime2 = TimeStamp(1) / 10^6;
time_whisk = (0:length(whisk)-1) * dt_whisk + starttime2;
timediff = starttime2 - starttime1;
% sniff = linterp(time_sniff,sniff,time_whisk);

% Resample
[p q] = rat(dsr/sr_sniff);    % resample respiration at 1000 Hz
sniff = resample(sniff,p,q);
[p q] = rat(dsr/sr_whisk);    % resample EMG at 1000 Hz
whisk = resample(whisk,p,q);
sr = dsr;  % new sampling rate
if timediff > 0.0001       % deal with misaligned data
    if ~isequal(sr,1000)
        error('Data realignment only works for samplin rate = 1000 Hz.')
    end
    leadtime = timediff * sr;    % lead of sniff in ms
    sniff(1:round(leadtime)) = [];  % drop lead: now they co-start
elseif timediff < -0.0001
    if ~isequal(sr,1000)
        error('Data realignment only works for samplin rate = 1000 Hz.')
    end
    lagtime = (-1) * timediff * sr;    % lag of sniff in ms
    whisk(1:round(lagtime)) = [];  % drop lead of whisk: now they co-start
end
ln = min(length(sniff),length(whisk));
sniff = sniff(1:ln);
whisk = whisk(1:ln);   % now they co-terminate

% Deal with individual artifacts
if isequal(rdate,'2008-09-03_16-42-14')
    whisk = whisk(1:2.25*10^5);
    sniff = sniff(1:2.25*10^5);
elseif isequal(rdate,'2008-09-21_21-08-00')
    whisk = whisk(1:1.54*10^6);
    sniff = sniff(1:1.54*10^6);
end

% Filter
nqf = sr / 2;
flt = fir1(1024,30/nqf,'low');
fsniff = filtfilt(flt,1,sniff);     % lowpass filter sniffing at 30 Hz
flt = fir1(4096,10/nqf,'high');
fwhisk = filtfilt(flt,1,whisk);     % highpass filter whisking at 10 Hz

% Respiration phase
resp_phase = angle(hilbert(standardize(fsniff)));    % Hilbert-transform

% Downsampled resp. phase
lenw = length(fwhisk);
lle = floor(lenw/wl) * wl;
resp_phase2 = resp_phase(1:lle);
resp_phase2 = resp_phase2(round(wl/2):wl:end);   % downsampled resp. phase aligned to whisking RMS

% Detect whisking
[st nd pk0 pk fwhisk_ rms_whisk thr] = sspw3(whisk,sr,wl,10,inf,[],[],thr);  % whisking detection, see SSPW3
if plotting
    figure
    plot(standardize(fwhisk))
    hold on
    line([pk; pk],[zeros(size(pk))-10; ones(size(pk))+20],'Color','green')
end

% Get resp. cycles
fn0 = valuecrossing(1:length(resp_phase),resp_phase',0,'down');   % limits for resp. cycles
fn20 = valuecrossing(1:length(resp_phase),resp_phase',0,'up');   % end of inhalation
fn2 = nan(size(fn0));
for k = 1:length(fn0)-1
    vcs = fn20(fn20>fn0(k)&fn20<fn0(k+1));
    fn2(k) = round(mean(vcs)/wl);
end
fn = round(fn0/wl);

if plotting
    figure
    line([fn*wl; fn*wl],[ones(size(fn))-4; ones(size(fn))+3],'Color','green')
    line([fn2*wl; fn2*wl],[ones(size(fn2))-4; ones(size(fn2))+3],'Color','red')
    hold on
    plot(resp_phase)
    
    figure
    line([fn*wl; fn*wl],[ones(size(fn))-4; ones(size(fn))+3],'Color','green')
    line([fn2*wl; fn2*wl],[ones(size(fn2))-4; ones(size(fn2))+3],'Color','red')
    hold on
    plot(standardize(fsniff))
end

% Sort whisking events
vdisc = round(pk/wl);
grp = zeros(size(vdisc));     % codes for types of whisking events 
time_grp = zeros(size(resp_phase2));   % codes time-wise 
time_grp2 = zeros(size(resp_phase2));   % indicator of cycles (time) before 1:2 whisking events
for k = 1:length(vdisc)
    ffv = find(fn<=vdisc(k),1,'last');    % fn: limits of respiration cycles
    if isempty(ffv) || ffv < 4 || ffv > length(fn)-4
        continue
    end
    fn_4pre = fn(ffv-3);     % cycle limit 4 before current events 
    fn_pppre = fn(ffv-2);    % cycle limit 3 before current events
    fn_ppre = fn(ffv-1);     % cycle limit 2 before current events
    fn_pre = fn(ffv);        % cycle limit 1 before current events
    fn_post = fn(ffv+1);     % cycle limit 1 after current events
    fn_ppost = fn(ffv+2);    % cycle limit 2 after current events
    fn_pppost = fn(ffv+3);   % cycle limit 3 after current events
    fn_4post = fn(ffv+4);    % cycle limit 4 after current events
    lv_pppre = length(vdisc(vdisc>=fn_4pre&vdisc<fn_pppre));    % whisking events in 3 cycles before
    lv_ppre = length(vdisc(vdisc>=fn_pppre&vdisc<fn_ppre));     % whisking events in 2 cycles before
    lv_pre = length(vdisc(vdisc>=fn_ppre&vdisc<fn_pre));        % whisking events in 1 cycle before
    lv = length(vdisc(vdisc>=fn_pre&vdisc<fn_post));            % whisking events in current cycle
    lv_post = length(vdisc(vdisc>=fn_post&vdisc<fn_ppost));     % whisking events in 1 cycle after
    lv_ppost = length(vdisc(vdisc>=fn_ppost&vdisc<fn_pppost));  % whisking events in 2 cycles after
    lv_pppost = length(vdisc(vdisc>=fn_pppost&vdisc<fn_4post)); % whisking events in 3 cycles after
    switch lv
        case 2
            grp(k) = 2;     % 2:1 coupling
            time_grp(floor(fn_pre):ceil(fn_post)) = 2;
        case 3              % 3:1 coupling
            grp(k) = 3;
            time_grp(floor(fn_pre):ceil(fn_post)) = 3;
        case 1              % 1:1 coupling
            if lv_pre >= 1 || lv_post >= 1   % at least one event in 1 cycle pre or post
                grp(k) = (1);
                time_grp(floor(fn_pre):ceil(fn_post)) = 1;
            else            % 1:2 coupling
                if lv_ppre >= 1 || lv_ppost >= 1    % no event in pre and post, but at least one event 2 cycles away
                    grp(k) = 0.5;
                    time_grp(floor(fn_ppre):ceil(fn_ppost)) = 0.5;
                    time_grp2(floor(fn_ppre):ceil(fn_pre-1)) = 1;   % indicator for cyles before and after the ones with the whisking events in 1:2 coupling
                    time_grp2(floor(fn_post):ceil(fn_ppost-1)) = 1;
                    if lv_ppre == 0     % prevent concatenation of 'out-of-phase' segments
                        time_grp(floor(fn_ppre):ceil(fn_ppre)) = 0;
                    elseif lv_ppost == 0
                        time_grp(floor(fn_ppost):ceil(fn_ppost)) = 0;
                    end
                else        % 1:3 coupling
                    if lv_pppre >= 1 || lv_pppost >= 1   % no event in 2 cycles around, but at least one event 3 cycles away
                        grp(k) = 1/3;
                        time_grp(floor(fn_pppre):ceil(fn_pppost)) = 1/3;
                        if lv_pppre == 0     % prevent concatenation of 'out-of-phase' segments
                            time_grp(floor(fn_pppre):ceil(fn_pppre)) = 0;
                        elseif lv_pppost == 0
                            time_grp(floor(fn_pppost):ceil(fn_pppost)) = 0;
                        end
                    end
                end
            end
    end
end

% Eliminate short segments
dfs = find(diff(time_grp));
segments = struct([]);
for k = 1:length(dfs)-1     % loop through uniform segments
    inx1 = dfs(k) + 1;  % segment start
    inx2 = dfs(k+1);    % segment end
    cseg = time_grp(inx1:inx2);   % segment (time domain)
    tcseg = length(cseg) / sr * wl;   % length of the segment in seconds
    cgrp = grp(vdisc>=inx1&vdisc<=inx2&grp~=0);   % whisking events within the segment
    if tcseg < 0.5 || length(cgrp) < 3  % exculde segments < 0.5 s and segments with < 3 whisking events
        time_grp(inx1:inx2) = 0;
        time_grp2(inx1:inx2) = 0;
        grp(vdisc>=inx1&vdisc<=inx2) = 0;
    else
        whs = find(vdisc>=inx1&vdisc<=inx2);
        vhs = vdisc(whs);
        whs2 = find(vdisc_rawwhisk>=inx1&vdisc_rawwhisk<=inx2);
        vhs2 = vdisc_rawwhisk(whs2);
        inst = fn(fn>=inx1&fn<=inx2+1);
        innd = nan(1,length(inst)-1);
        for ks = 1:length(inst)-1
            ksinx = find(fn2>=inst(ks),1,'first');
            innd(ks) = fn2(ksinx);
            if isequal(inst(ks),innd(ks))
                if fn2(ksinx+1) < inst(ks+1)
                    innd(ks) = fn2(ksinx+1);
                end
            end
        end
        exst = innd;
        exnd = inst(2:end);
        inst = inst(1:end-1);
        frq = length(exst) / tcseg;
        if ~b_isconstant(grp(whs(1:end-1)))
            error('Technical error 236.')
        end
        if ~isequal(length(inst),length(innd),length(exst),length(exnd))
            error('Technical error 244.')
        end
        cgrp = grp(whs(1));
        seginx = length(segments) + 1;
        segments(seginx).rat = rrat;
        segments(seginx).session = rdate;
        segments(seginx).start = inx1;
        segments(seginx).end = inx2;
        segments(seginx).length = tcseg;
        segments(seginx).frequency = frq;
        segments(seginx).group = cgrp;
        segments(seginx).inhalation_start = inst;
        segments(seginx).inhalation_end = innd;
        segments(seginx).exhalation_start = exst;
        segments(seginx).exhalation_end = exnd;
        segments(seginx).whisk_events = vhs;
        segments(seginx).raw_whisk_events = vhs2;
    end
end

% Split 1:1 coupling
grp_split = grp;   % mark every odd whisk in 1:1 coupling bouts
time_grp_split = time_grp;   % same in time
for k = 1:length(grp)-1
    ffv = find(fn<vdisc(k),1,'last');
    fn_pre = fn(ffv);        % cycle limit 1 before current events
    fn_post = fn(ffv+1);     % cycle limit 1 after current events
    if grp(k) == 1 && (grp(k-1) ~= 1 || grp_split(k-1) == 1)
        grp_split(k) = 0.9;   % code for odd whisk events in 1:1 coupling bouts
        time_grp_split(floor(fn_pre):ceil(fn_post)) = 0.9;
    end
end

% Respiration phase - whisking relationship
time_grp_orig = time_grp;
rms_whisk_orig = rms_whisk;
rms_whisk(rms_whisk<0) = 0;     % rectify RMS: due to filtering, below-0 values could appear
time_grp1 = time_grp;   % 1:1 coupling
time_grp1(time_grp~=1) = 0;     % indicator variable for 1:1 coupling (time)
time_grp1(time_grp==1) = 1;
[R p] = dphist(time_grp1,rms_whisk,resp_phase2);    % whisking RMS conditioned on resp. phase
fnm = [resdir 'PHASE1TO1_' rrat '_' rdate '.fig'];   % save figure
saveas(gcf,fnm)
fnm = [resdir 'PHASE1TO1_' rrat '_' rdate '.mat'];   % save variables
save(fnm,'R','p','time_grp1')

time_grp1 = time_grp;   % 1:2 coupling
time_grp1(time_grp~=0.5) = 0;     % indicator variable for 1:2 coupling (time)
time_grp1(time_grp==0.5) = 1;
time_grp1(time_grp2==1) = 0;    % disregard those not having whisking
[R p] = dphist2(time_grp1,time_grp2,rms_whisk,resp_phase2);    % whisking RMS conditioned on resp. phase
fnm = [resdir 'PHASE1TO2_' rrat '_' rdate '.fig'];   % save figure
saveas(gcf,fnm)     % plot shows every cycle before/after the detected events in black, see 'time_grp2'
fnm = [resdir 'PHASE1TO2_' rrat '_' rdate '.mat'];   % save variables
save(fnm,'R','p','time_grp1','time_grp2')

time_grp1 = time_grp;   % 2:1 coupling
time_grp1(time_grp~=2) = 0;     % indicator variable for 2:1 coupling (time)
time_grp1(time_grp==2) = 1;
[R p] = dphist(time_grp1,rms_whisk,resp_phase2);    % whisking RMS conditioned on resp. phase
fnm = [resdir 'PHASE2TO1_' rrat '_' rdate '.fig'];   % save figure
saveas(gcf,fnm)     % plot shows every cycle before the detected events in black, see 'time_grp2'
fnm = [resdir 'PHASE2TO1_' rrat '_' rdate '.mat'];   % save variables
save(fnm,'R','p','time_grp1')

% 1:1 couling, split into two groups to compare with 1:2
time_grp_split1 = time_grp_split;   % 1:1 coupling, split
time_grp_split2 = time_grp_split;
time_grp_split1(time_grp_split~=0.9) = 0;     % indicator variable for odd cycles (time)
time_grp_split1(time_grp_split==0.9) = 1;
time_grp_split2(time_grp_split~=1) = 0;     % indicator variable for even cycles (time)
time_grp_split2(time_grp_split==1) = 1;
[R p] = dphist2(time_grp_split1,time_grp_split2,rms_whisk,resp_phase2);    % whisking RMS conditioned on resp. phase
fnm = [resdir 'PHASE1TO1SPLIT_' rrat '_' rdate '.fig'];   % save figure
saveas(gcf,fnm)     % plot shows every even cycle in black
fnm = [resdir 'PHASE1TO1SPLIT_' rrat '_' rdate '.mat'];   % save variables
save(fnm,'R','p','time_grp_split1','time_grp_split2')

% Save additional variables
fnm = [resdir 'PHASE_' rrat '_' rdate '.mat'];   % save variables
save(fnm,'rms_whisk','rms_whisk_orig','vdisc','thr','resp_phase2','fn','time_grp','grp','segments',...
    'wl','csc1','csc2','msgn','rrat','rdate')

% keyboard

% -------------------------------------------------------------------------
function [st nd pk1 pk2 feeg rms thr] = sspw3(eeg,sr,wl,filtL,filtU,mlt1,mlt2,thr)
%SSPW3   Wave detection by RMS thresholding.
%   [ST ND PK1 PK2 FEEG RMS THR] = SSPW3(EEG,SR,WL,L,U,M1,M2) detectects
%   waves (filtered between L and U Hz) in EEG sampled at SR and returns
%   starting, end and peak points of epochs in ST, ND, PK1 and PK2. PK1
%   localizes RMS peaks while PK2 finds the peaks of the original signal.
%   SSPW2 uses the following criteria for wave detection: root-mean-square
%   (window size determined by WL) of filtered EEG should reach mean(RMS) +
%   M2 * std(RMS) and peak RMS should reach mean(RMS) + M1 * std(RMS). If
%   M1 is empty, the user is prompted to assess a threshold manually. If M2
%   is empty, 0.1 of the discrimination thresold is used. Filtered EEG
%   (FEEG), root-mean-square (RMS) and discrimination threshold (THR) are
%   also returned.
%
%   SSPW3(EEG,SR,WL,L,U,M1,M2,THR) accepts a predifined threshold for
%   assessing RMS peaks. THR overwrites M1.
%
%   In SSPW3, RMS is high-pass filtered at 50 Hz.
%
%   See also SSPW.

% Input argument check
error(nargchk(7,8,nargin))
if nargin < 8
    thr = [];
end

% Filtering
nqf = sr / 2;
if filtU == inf
    b = fir1(2048,filtL/nqf,'high');
elseif filtL == 0
    b = fir1(2048,filtU/nqf,'low');
else
    b = fir1(2048,[filtL filtU]/nqf);
end
feeg = filtfilt(b,1,eeg);

% Root mean square
leneeg = length(feeg);
% wl = 100;
lle = floor(leneeg/wl) * wl;
feeg2 = reshape(feeg(1:lle),wl,lle/wl);
prms = sqrt(sum(feeg2.^2)) / sqrt(wl);
if length(prms) < 3 * 2048
    warning('sspw3:filterOrder','Lower filter order.')
    fo = 2 ^ round(log2(length(prms)/3));
    disp(['Filter order: ' num2str(fo)])
else
    fo = 2048;
end
b = fir1(fo,50/nqf,'high');       % highpass filter
rms = filtfilt(b,1,prms);

% Discriminate RMS: RMS peak during sharpwave should reach mean(RMS) + mlt * std(RMS)
mrms = mean(rms);
sdrms = std(rms);
if ~isempty(mlt1)
    if isempty(thr)
        thr = mrms + mlt1 * sdrms;
    end
    pks = disc(rms,thr);
else
    if isempty(thr)
        [pks thr] = b_udisc(rms);
    else
        pks = disc(rms,thr);
    end
end

% Ripple start, end: RMS should cross mean(RMS) + std(RMS)
lenp = length(pks);
st = zeros(1,lenp);
nd = zeros(1,lenp);
pk1 = zeros(1,lenp);
pk2 = zeros(1,lenp);
mrms = mean(rms);
sdrms = std(rms);
if ~isempty(mlt2)
    v = mrms + mlt2 * sdrms;
else
    v = thr * 0.1;
end
lup = rms < v & [rms(2:end) v] > v;   % value-crossings
lup2 = find(lup) * wl;
ldown = rms < v & [v rms(1:end-1)] > v;
ldown2 = find(ldown) * wl;
pks = pks * wl;
for k = 1:length(pks)
    sst = lup2(find(lup2<pks(k),1,'last'));
    nnd = ldown2(find(ldown2>pks(k),1,'first'));
    if isempty(sst) || isempty(nnd)
        continue
    end
    st(k) = sst;
    nd(k) = nnd;
    pk1(k) = sst + (find(rms(sst/wl:nnd/wl)==max(rms(sst/wl:nnd/wl)))-1) * wl;
    pk2(k) = sst + find(feeg(sst:nnd)==max(feeg(sst:nnd)));
end
st = st(st>0);
nd = nd(nd>0);
[st m1] = unique(st);
[nd m2] = unique(nd);
[pk1 m3] = unique(pk1);
[pk2 m4] = unique(pk2);
if ~isequal(m1,m2,m3,m4)
    error('Technical error 66.')
end

% -------------------------------------------------------------------------
function [R p] = dphist(time_grp1,rms_whisk,resp_phase2) 
edges = -pi:2*pi/18:pi;     % phase histogram bin limits
cnts = (edges(1:end-1) + edges(2:end)) / 2;     % phase histogram bin centers
sfwhisk = rms_whisk' .* time_grp1;
zeroinx = time_grp1==0;
spvr = [resp_phase2 sfwhisk];
spvr(zeroinx,:) = [];
spvr = sortrows(spvr,1);
[mn_phase mn_wh] = phasedep(spvr,edges);
figure
plot([cnts cnts+2*pi],[mn_wh mn_wh],'r')      % conditional distribution: E(whisk|phi1<phase<phi2)
if ~isempty(find(time_grp1,1))
    [R p] = lincirc_corr2(mn_wh',cnts')
else
    R = NaN;
    p = NaN;
end

% -------------------------------------------------------------------------
function [R p] = dphist2(time_grp1,time_grp2,rms_whisk,resp_phase2) 
edges = -pi:2*pi/18:pi;     % phase histogram bin limits
cnts = (edges(1:end-1) + edges(2:end)) / 2;     % phase histogram bin centers
sfwhisk = rms_whisk' .* time_grp1;
zeroinx = time_grp1==0;
spvr = [resp_phase2 sfwhisk];
spvr(zeroinx,:) = [];
spvr = sortrows(spvr,1);
[mn_phase mn_wh] = phasedep(spvr,edges);
figure
plot([cnts cnts+2*pi],[mn_wh mn_wh],'r')      % conditional distribution: E(whisk|phi1<phase<phi2)
if ~isempty(find(time_grp1,1))
    [R p] = lincirc_corr2(mn_wh',cnts')
else
    R = NaN;
    p = NaN;
end

sfwhisk = rms_whisk' .* time_grp2;
zeroinx = time_grp2==0;
spvr = [resp_phase2 sfwhisk];
spvr(zeroinx,:) = [];
spvr = sortrows(spvr,1);
[mn_phase mn_wh] = phasedep(spvr,edges);
hold on
plot([cnts cnts+2*pi],[mn_wh mn_wh],'k')     % overlay a second phase histogram