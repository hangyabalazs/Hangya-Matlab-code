function [st nd feeg] = sspw(eeg,sr)
%SSPW   Sharpwave detection.
%   [ST ND] = SSPW(EEG,SR) detectects ripples in EEG sampled at SR and
%   returns starting and end points of ripples in ST and ND, respectively.
%   SSPW uses the following criteria for ripples: root-mean-square of EEG 
%   during a ripple should reach mean(RMS) + std(RMS) and peak RMS should
%   reach mean(RMS) + 7 * std(RMS).
%
%   See also SSPW_NONTHETA.

% Input argument check
error(nargchk(2,2,nargin))

% Filtering
nqf = sr / 2;
b = fir1(2048,[90 140]/nqf);
feeg = filtfilt(b,1,eeg);

% Root mean square
leneeg = length(feeg);
wl1 = 100;
lle1 = floor(leneeg/wl1) * wl1;
feeg21 = reshape(feeg(1:lle1),wl1,lle1/wl1);
rms1 = sqrt(sum(feeg21.^2)) / sqrt(wl1);
wl2 = 100;
lle2 = floor(leneeg/wl2) * wl2;
feeg22 = reshape(feeg(1:lle2),wl2,lle2/wl2);
rms2 = sqrt(sum(feeg22.^2)) / sqrt(wl2);

% Discriminate RMS: RMS peak during sharpwave should reach mean(RMS) + 5 * std(RMS)
mrms1 = mean(rms1);
sdrms1 = std(rms1);
thr = mrms1 + 5 * sdrms1;
pks = disc(rms1,thr);

% Ripple start, end: RMS should cross mean(RMS) + std(RMS)
lenp = length(pks);
st = zeros(1,lenp);
nd = zeros(1,lenp);
mrms2 = mean(rms2);
sdrms2 = std(rms2);
v = mrms2 + sdrms2;
lup = rms2 < v & [rms2(2:end) v] > v;   % value-crossings
lup2 = find(lup) * wl2;
ldown = rms2 < v & [v rms2(1:end-1)] > v;
ldown2 = find(ldown) * wl2;
pks = pks * wl1;
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