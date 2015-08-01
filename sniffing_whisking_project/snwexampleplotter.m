function snwexampleplotter(cseg)
%SNWEXAMPLEPLOTTER   Plot example for sniffing-whisking project.
%   SNWEXAMPLEPLOTTER(CSEG) plots sniffing (raw, filtered, phase) and
%   whisking (high-pass filtered, RMS, detected whisk events) for a segment
%   (SEG, see SNWSNIFFPLOT2 for details on the 'segment' structure).
%
%   See also SNWSNIFFPLOT2 and SNWEXAMPLEPLOTTER2.

% Import session data from Excel
global DATAPATH
tblfile = [DATAPATH 'SniffWhisk\data_selection.xlsx'];
[tbl0 tbl] = xlsread(tblfile);
global DATADIR
ptblinx = find(strcmp(tbl(:,1),cseg.rat));
ptblinx2 = find(strcmp(tbl(ptblinx,2),cseg.session));
tblinx = ptblinx(ptblinx2);

% Data import
plotting = false;
dsr = 1000;   % downsample rate
rrat = cseg.rat;
rdate = cseg.session;
csc1 = tbl{tblinx,3};
csc2 = tbl{tblinx,4};
msgn = tbl0(tblinx,1);
thr = tbl0(tblinx,2);
wl = tbl0(tblinx,3);

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

% Root mean square of whisking
lenw = length(fwhisk);
lle = floor(lenw/wl) * wl;
wh2 = reshape(fwhisk(1:lle),wl,lle/wl);
urms_whisk = sqrt(sum(wh2.^2)) / sqrt(wl);
urms_whisk = urms_whisk';

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

% Whisk events
vdisc = round(pk/wl);

% keyboard

% plot
segstinx = round((cseg.start-0.2)*sr);
segndinx = round((cseg.end+0.2)*sr);
lfsniff = fsniff(segstinx:segndinx);
lsniff = sniff(segstinx:segndinx);
lfwhisk = fwhisk(segstinx:segndinx);
lwhisk = whisk(segstinx:segndinx);
lrms_whisk = urms_whisk(round(segstinx/wl):round(segndinx/wl));
lresp_phase = resp_phase(segstinx:segndinx);
lfn2 = fn2(fn2*wl>segstinx&fn2*wl<segndinx) * wl - segstinx;
lvdisc = vdisc(vdisc*wl>segstinx&vdisc*wl<segndinx) * wl - segstinx;
figure
plot(standardize(lfsniff),'k')
hold on
plot(standardize(lsniff),'r')
hold on
plot(lresp_phase,'Color',[0.7 0.7 0.7])
line([lfn2; lfn2],[zeros(size(lfn2))-pi; zeros(size(lfn2))+pi],'Color',[0 0.7 0.3])
% set(gcf,'Position',[1924 630 1900 342])
setappdata(gcf,'segment',cseg)
S1 = gca;

figure
plot(lfwhisk/300,'r')
hold on
plot(linspace(1,length(lfwhisk),length(lrms_whisk)),lrms_whisk/300,'k')
line([lvdisc; lvdisc],[zeros(size(lvdisc))-1; zeros(size(lvdisc))+50],'Color',[0 0.7 0.3])
% set(gcf,'Position',[1924 71 1900 511])
setappdata(gcf,'segment',cseg)
S2 = gca;

keyboard

linkaxes([S1 S2],'x')



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