function [st nd pk feeg rms thr] = sspw3(eeg,sr,wl,filtL,filtU,mlt1,mlt2)
%SSPW3   Wave detection by RMS thresholding.
%   [ST ND  PK FEEG RMS THR] = SSPW3(EEG,SR,WL,L,U,M1,M2) detectects waves
%   (filtered between L and U Hz) in EEG sampled at SR and returns
%   starting, end and peak points of epochs in ST, ND and PK, respectively.
%   SSPW2 uses the following criteria for wave detection: root-mean-square
%   (window size determined by WL) of filtered EEG should reach mean(RMS) +
%   M2 * std(RMS) and peak RMS should reach mean(RMS) + M1 * std(RMS). If
%   M1 is empty, the user is prompted to assess a threshold manually. If M2
%   is empty, 0.1 of the discrimination thresold is used. Filtered EEG
%   (FEEG), root-mean-square (RMS) and discrimination threshold (THR) are
%   also returned.
%
%   In SSPW3, RMS is high-pass filtered at 50 Hz.
%
%   See also SSPW.

% Input argument check
error(nargchk(7,7,nargin))

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
    thr = mrms + mlt1 * sdrms;
    pks = disc(rms,thr);
else
    [pks thr] = b_udisc(rms);
end

% Ripple start, end: RMS should cross mean(RMS) + std(RMS)
lenp = length(pks);
st = zeros(1,lenp);
nd = zeros(1,lenp);
pk = zeros(1,lenp);
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
    pk(k) = sst + (find(rms(sst/wl:nnd/wl)==max(rms(sst/wl:nnd/wl)))-1) * wl;
end
st = st(st>0);
nd = nd(nd>0);
[st m1] = unique(st);
[nd m2] = unique(nd);
[pk m3] = unique(pk);
if ~isequal(m1,m2,m3)
    error('Technical error 66.')
end