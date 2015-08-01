function [pow,f,newstep] = b_powercall_for_applications
%POWERCALL_FOR_APPLICATIONS   Calls wavelet for EEG.
%   POWERCALL_FOR_APPLICATIONS calls wavelet for EEG with fixed parameters (e.g.
%   newstep for downsampling is 25, which means 400 Hz downsampling using 10 kHz 
%   sampled data). Makes calling wavelet easier from functions.
%
%   Unlike WAVELETCALL_FOR_APPLICATION, computes only power.
%
%   POWERCALL_FOR_APPLICATIONS calls WAVELET_POW3.
%
%   See also WAVELET, WAVELET_POW3 and WAVELETCALL_FOR_APPLICATIONS.

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
[pow,period,scale,coi] = b_wavelet_pow3(sst,dt,pad,dj,s0,j1,mother,param,mis);