function [pow,phase,f] = eegwavelet(dat,sr)
%EEGWAVELET   Wavelet calculation.
%   [POW,PHASE,F] = EEGWAVELET(DATA,SR) performs wavelet calculation on
%   input EEG data (DATA) sampled at SR sampling rate. Wavelet power (POW),
%   phase (PHASE) and scale vector (F) are returned. DATA is first
%   standardized. A minimal interesting frequency of 0.5 Hz is set.
%
%   See also UNITWAVELET.

% Prepare for wavelet transformation
variance = std(dat) ^ 2;
dat = (dat - mean(dat)) / sqrt(variance) ;
n = length(dat);
dt = 1 / sr;
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
[wave,period,scale,coi] = b_wavelet_new3(dat,dt,pad,dj,s0,j1,mother,param,mis);
pow = abs(wave).^2;
phase = angle(wave);