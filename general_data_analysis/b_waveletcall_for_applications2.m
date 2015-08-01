function [wave_eeg,wave_unit,f,newstep] = b_waveletcall_for_applications2
%WAVELETCALL_FOR_APPLICATIONS2   Calls wavelet for EEG and unit.
%   WAVELETCALL_FOR_APPLICATIONS calls wavelet for both EEG and unit with fixed 
%   parameters (e.g. newstep for downsampling is 25, which means 400 Hz downsampling
%   using 10 kHz sampled data). Unit is first convolved with sinc function (see 
%   SINCCONV2 for details). Makes calling wavelet easier from functions.
%
%   WAVELETCALL_FOR_APPLICATIONS calls WAVELET_NEW3.
%
%   See also WAVELET, WAVELET_NEW2, WAVELET_NEW3, WAVELETCALL_FOR_APPLICATIONS
%   and SINCCONV2.

% Create variables in workspace
b_var2ws('all','caller')

% Create variables for wavelet
% wholeeeg = data(datinx1:datinx2,1);
newstep = 25;
resamp = 10000/newstep;

sst_eeg = eeg(1:newstep:end);   % downsample eeg
sst_unit = b_sincconv2(unit,vdisc,newstep);    % unit sinc convolution

% Standardization
sst_eeg = (sst_eeg - mean(sst_eeg)) / std(sst_eeg);
sst_unit = (sst_unit - mean(sst_unit)) / std(sst_unit);

% Wavelet transformation
dt = 1 / resamp;    %resample on 400 Hz
n = length(sst_eeg);
time = [0:n-1] * dt + 0;
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
[wave_eeg,period,scale,coi] = b_wavelet_new3(sst_eeg,dt,pad,dj,s0,j1,mother,param,mis);
[wave_unit,period,scale,coi] = b_wavelet_new3(sst_unit,dt,pad,dj,s0,j1,mother,param,mis);