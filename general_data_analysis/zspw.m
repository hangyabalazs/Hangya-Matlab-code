function [st nd mx feeg] = zspw(eeg,sr,thmp)
%ZSPW   Sharpwave detection.
%   [ST ND] = ZSPW(EEG,SR) detectects ripples in EEG sampled at SR and
%   returns starting and end points of ripples in ST and ND, respectively.
%   ZSPW uses the following criteria for ripples: root-mean-square of EEG 
%   during a ripple should reach mean(RMS) + std(RMS) and peak RMS should
%   reach mean(RMS) + 5 * std(RMS).
%   ZSPW plots a figure with the original EEG, the filtered EEG, the RMS
%   and the applied thresholds. It is possible to zoom on the ripples using
%   the 'x', 'y' and 'o' buttons (see ZSPW_KEYPRESS).
%
%   [ST ND] = ZSPW(EEG,SR,M) accepts an input argument (M) determining the
%   multiplier of the RMS SD for thresholding.
%
%   [ST ND MX] = ZSPW(EEG,SR,M) returns ripple peaks in MX. Peaks are
%   defined as maximum locations of the filtered (90-200 Hz) EEG.
%
%   [ST ND MX FEEG] = ZSPW(EEG,SR) returns EEG band-pass filtered between 90
%   and 200 Hz in FEEG.
%
%   See also SSPW and ZSPW_KEYPRESS.

% Input argument check
error(nargchk(2,3,nargin))
if nargin < 3
    thmp = 5;
end

% Filtering
nqf = sr / 2;
b = fir1(2048,[90 200]/nqf);
feeg = filtfilt(b,1,eeg);

% Root mean square
leneeg = length(feeg);
wl = sr / 100;     % 10 ms window
lle = floor(leneeg/wl) * wl;
feeg2 = reshape(feeg(1:lle),wl,lle/wl);
rms = sqrt(sum(feeg2.^2)) / sqrt(wl);

% Discriminate RMS: RMS peak during sharpwave should reach mean(RMS) + 5 * std(RMS)
mrms = mean(rms);
sdrms = std(rms);
thr = mrms + thmp * sdrms;
pks = disc(rms,thr);

% Ripple start, end: RMS should cross mean(RMS) + std(RMS)
lenp = length(pks);
st = zeros(1,lenp);
nd = zeros(1,lenp);
v = mrms + sdrms;
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
end
st = st(st>0);
nd = nd(nd>0);
[st m1] = unique(st);
[nd m2] = unique(nd);
if ~isequal(m1,m2)
    error('Technical error 50.')
end
lenst = length(st);
mx = zeros(1,lenst);    % ripple peaks: maximum location of the filtered EEG
for k = 1:lenst
    mxv = max(feeg(st(k):nd(k)));
    mxl = find(feeg(st(k):nd(k))==mxv);
    mx(k) = st(k) + mxl(1);
end

% Plot
times = linspace(0,length(eeg)/sr,length(eeg));
rms_times = linspace(0,length(eeg)/sr,length(rms));
H = figure;
plot(times,eeg/3)
hold on
plot(times,feeg,'c')
plot(rms_times,rms,'r')
line([times(1) times(end)],[thr thr],'Color','g')
line([times(1) times(end)],[v v],'Color','g','LineStyle',':')
for k = 1:lenst
    plot(times(st(k):nd(k)),eeg(st(k):nd(k))/3,'m')
end
kpf = 'zspw_keypress';   % set keypress function
set(H,'KeyPressFcn',kpf)
x_lim = ylim;   % set application data
setappdata(H,'x_lim',[times(1) times(end)])
setappdata(H,'y_lim',ylim)
rippleindex = 0;
setappdata(H,'rippleindex',rippleindex)
setappdata(H,'st',st)
setappdata(H,'nd',nd)
setappdata(H,'sr',sr)