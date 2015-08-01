function [rWCo cWCo] = icoh(eeg1,eeg2,sr)
%ICOH   Calculates wavelet coherence.
%   [RW CW] = ICOH(EEG1,EEG2,SR) calculates wavelet coherence (RW) between
%   EEG1 and EEG2 (both sampled on SR). Control data (CW) is calculated
%   between EEG1 and segment-shuffled EEG2.
%
%   Selected frequency band: 0.5 - 4 Hz.
%
%   Reference: Lachaux JP, Lutz A, Rudrauf D, Cosmelli D, Queyen MLV,
%   Martinerie J, Varela F (2002) Estimating the time-course of coherence
%   betweeen single-trial brain signals: an introduction to wavelet
%   coherence. Neurophysiol Clin 32:157-174.
%
%   See also WAVELET_NEW3.

% Input argument check
error(nargchk(3,3,nargin))

% Create random eeg (segment shuffle)
seglen = length(eeg2) * sr / 1000;
ic1 = fix(seglen/sr);
increm = round(seglen/ic1);     % 1-1.2 sec. long segments
ed = {};
for t = 1:ic1
    ind1 = (t - 1) * increm + 1;
    ind2 = ind1 + increm - 1;
    ed{end+1} = eeg2(ind1:ind2);
end
led = length(ed);
rp = randperm(led);
while any(rp==[1:led])
    rp = randperm(led);
end
psed = [];
for j = 1:led
    psed = [psed ed{rp(j)}];
end

% Wavelet transformation of eeg
wavea = eeg_wavelet(eeg1);    % eeg1
waveb = eeg_wavelet(eeg2);   % eeg2
[wavec ,f] = eeg_wavelet(psed);    % random eeg

% Wavelet coherence
fnd = find(f>4);    % frequency band bounderies
pwind1 = fnd(end);
wlen = min(length(psed),seglen);
wave1 = wavea(pwind1:end,1:wlen);  % eeg1
wave2 = waveb(pwind1:end,1:wlen);  % eeg2
wave3 = wavec(pwind1:end,1:wlen);  % random eeg
f = f(pwind1:end);

rWCo = wavelet_coherence(wave1,wave2,f,sr);    % real
cWCo = wavelet_coherence(wave1,wave3,f,sr);    % control
clear wavea waveb wavec



% -------------------------------------------------------------------------
% WAVELET
% -------------------------------------------------------------------------
function [wave,f] = eeg_wavelet(dat)

% Prepare for wavelet transformation
variance = std(dat) ^ 2;
dat = (dat - mean(dat)) / sqrt(variance) ;
n = length(dat);
dt = 1 / 1000;
pad = 1;
dj = 0.08;    
j1 = ceil((1/dj) * log2(n/2));
j1 = ceil(j1);
j = (0:j1);
s0 = 2 * dt; 
s = s0 .* 2 .^ (j * dj);
omega0 = 6;
c = 4 * pi / (omega0 + sqrt(2+omega0^2));
fperiod = c .* s;
f = 1 ./ fperiod;
lag1 = 0.72;
param = -1;
mif = 0.5;          %minimal interesting frequency
mis = find(f>mif);
mis = mis(end);     %maximal interesting scale
mother = 'Morlet';

% Wavelet transformation
[wave,period,scale,coi] = b_wavelet_new3(dat,dt,pad,dj,s0,j1,mother,param,mis);
f = f(f>mif);



% -------------------------------------------------------------------------
% COHERENCE
% -------------------------------------------------------------------------
function WCo = wavelet_coherence(Wx,Wy,f,samprate)

nsc = size(Wx,1);       % number of scales
[k1 k2] = size(Wx);
CWxx = Wx .* conj(Wx);
CWyy = Wy .* conj(Wy);
CWxy = Wx .* conj(Wy);

ncy = 10;       % >5 for reasonably small bias
SWxx = zeros(k1,k2);
SWyy = zeros(k1,k2);
SWxy = zeros(k1,k2);
for k = 1:nsc
    fr = f(k);
    dlt = ncy / fr;
    delta = dlt * samprate;
    delta = round(delta);
    delta = delta - rem(delta,2) + 1;   % making the integration interval length odd
    deltahalf = (delta - 1) / 2;
    flt = ones(1,delta);
    pS = conv(CWxx(k,:),flt);
    SWxx(k,:) = pS(deltahalf+1:end-deltahalf);
    pS = conv(CWyy(k,:),flt);
    SWyy(k,:) = pS(deltahalf+1:end-deltahalf);
    pS = conv(CWxy(k,:),flt);
    SWxy(k,:) = pS(deltahalf+1:end-deltahalf);
end
WCo = abs(SWxy) ./ (SWxx .* SWyy).^0.5;