function [pow,phase,f] = unitwavelet(vdisc,lenu,sr)
%UNITWAVELET   Wavelet calculation.
%   [POW,PHASE,F] = UNITWAVELET(DATA,LENU,SR) performs wavelet calculation 
%   on discriminated unit data (VDISC) sampled at SR sampling rate for a
%   data length of LENU. Wavelet power (POW), phase (PHASE) and scale
%   vector (F) are returned. DATA is first convolved by a sinc kernel and
%   standardized. A minimal interesting frequency of 0.5 Hz is set.
%
%   See also EEGWAVELET.

% Sinc convolution
fs = sr;     % unit
dto = 1 / fs;
ts = zeros(1,lenu);
ts(vdisc) = 1;
du = diff(vdisc);
fdu = 1 ./ du;
fdu = [fdu 1/sr];
fcut = 100; 
fsnew = 1000;
dtnew = 1 / 1000;
fsold = sr;
fsratio = fsnew / fsold;
told = vdisc * dto * fcut;
tnew = (1:lenu*fsratio) * dtnew * fcut;
lentold = length(told);
zint = 0;
for i = 1:lentold
    zint = zint + sinc(tnew-told(i));
end

% Prepare for wavelet transformation
variance = std(zint) ^ 2;
zint = (zint - mean(zint)) / sqrt(variance) ;
n = length(zint);
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
mif = 0.5;          %minimal intresting frequency
mis = find(f>mif);
mis = mis(end);     %maximal intristing scale
mother = 'Morlet';

% Wavelet transformation
[wave,period,scale,coi] = b_wavelet_new3(zint,dt,pad,dj,s0,j1,mother,param,mis);
pow = abs(wave).^2;
phase = angle(wave);