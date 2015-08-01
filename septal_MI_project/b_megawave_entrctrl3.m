function [wavef] = b_megawave_entrctrl3(eeg,vdisc)
%MEGAWAVE_ENTRCTRL3   Calculates wavelet transform for delayed unit.
%   EEG and 'vdisc' are needed as input arguments. The time series are preprocessed through sinc 
%   convolution before wavelet transformation.
%
%   F = MEGAWAVE_ENTRCTRL(EEG,VDISC) returns F output argument: delayed unit wavelet with 10 seconds
%   delay.
%
%   See also MEGAWAVE_ENTRCTRL and MEGAWAVERUN_ENTRCTRL3.

% Input arguments check
error(nargchk(2,2,nargin));

% Create random unit
lenu = length(eeg);
r = random('exp',1000,1,1500);
s = cumsum(r);
psvd = unique(ceil(s));     % 'pseudo vdisc'
psvd = psvd(find(psvd<lenu));

% Create delayed unit
flag = vdisc(find(vdisc>(lenu-100000)));
iflag = vdisc(find(vdisc<=(lenu-100000)));
dvd = [flag-lenu+100000 iflag+100000];     % 'delayed vdisc': last second of unit comes to first

% Sinc convolution & wavelet transformation
wavef = unit_wavelet(dvd,lenu);     % delayed unit wavelet



% --------------------------------------------------------------------------------------------
function wave = unit_wavelet(vdisc,lenu)

% Sinc convolution
fs = 10000;     % unit
dto = 1 / fs;
ts = zeros(1,lenu);
ts(vdisc) = 1;
du = diff(vdisc);
fdu = 1 ./ du;
fdu = [fdu 0.0001];
fcut = 100; 
fsnew = 1000;
dtnew = 1 / 1000;
fsold = 10000;
fsratio = fsnew / fsold;
told = vdisc * dto * fcut;
tnew = (1:lenu*fsratio) * dtnew * fcut;
lentold = length(told);
zint = 0;
for i = 1:lentold
    zint = zint + sinc(tnew-told(i));
end

% Prepare for wavelet transformation
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
[wave,period,scale,coi] = b_wavelet_pow3(zint,dt,pad,dj,s0,j1,mother,param,mis);