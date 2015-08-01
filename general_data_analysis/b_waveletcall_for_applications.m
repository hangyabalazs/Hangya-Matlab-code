function [wave,f,newstep] = b_waveletcall_for_applications
%WAVELETCALL_FOR_APPLICATIONS   Calls wavelet for EEG.
%   WAVELETCALL_FOR_APPLICATIONS calls wavelet for EEG with fixed parameters (e.g.
%   newstep for downsampling is 25, which means 400 Hz downsampling using 10 kHz 
%   sampled data). Makes calling wavelet easier from functions.
%
%   WAVELETCALL_FOR_APPLICATIONS calls WAVELET_NEW2.
%
%   See also WAVELET, WAVELET_NEW2, WAVELET_NEW3 and WAVELETCALL_FOR_APPLICATIONS2.

% Create variables in workspace
b_var2ws('in','caller')

% Create variables for wavelet
wholeeeg = data(datinx1:datinx2,1);
lenwe = length(wholeeeg);
newstep = 25;
resamp = 10000/newstep;
sst = wholeeeg(1:newstep:lenwe);

% Standardization
variance = std(sst)^2;
sst = (sst - mean(sst)) / sqrt(variance);

% Wavelet transformation
dt = 1 / resamp;    %resample on 400 Hz
time = [0:length(sst)-1] * dt + 0;
n = length(sst);
pad = 1;      
dj = 0.04;
s0 = 2 * dt;
j1 = ((1 / dj) * log2(n/2)) * 2;
j1 = ceil(j1);
j = (0:j1);
s = s0 .* 2 .^ (j * dj);
omega0 = 6;
c = 4 * pi / (omega0 + sqrt(2+omega0^2));
fperiod = c .* s;
f = 1 ./ fperiod;
lag1 = 0.72;  
mother = 'Morlet';
param = -1;
mif = 0.5;          %minimal intresting frequency
mis = find(f>mif);
mis = mis(end);     %maximal intristing scale
[wave,period,scale,coi] = b_wavelet_new2(sst,dt,pad,dj,s0,j1,mother,param,mis);