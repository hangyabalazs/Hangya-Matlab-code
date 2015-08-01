function [wavea,waveb,wavec,waved,wavee,wavef] = b_2megawave_entrctrl(eeg,vdisc)
%2MEGAWAVE_ENTRCTRL   Calculates wavelet transform for unit, eeg, random sample, random generated and delayed unit.
%   EEG and 'vdisc' are needed as input arguments. The time series are preprocessed through sinc 
%   convolution before wavelet transformation.
%
%   [A,B,C,D,E,F] = MEGAWAVE_ENTRCTRL(EEG,VDISC) returns 6 output arguments: eeg wavelet (A), unit 
%   wavelet (B), random sample wavelet (C), random generated unit wavelet (D), delayed unit wavelet
%   with 1 second delay (E), delayed unit wavelet with 10 seconds delay (F). Random generated unit
%   is a randomized Poisson process.
%
%   See also 2MEGAWAVERUN_ENTRCTRL.

% Input arguments check
error(nargchk(2,2,nargin));

% Create random unit
lenu = length(eeg);
r = random('exp',1000,1,1500);
s = cumsum(r);
psvd = unique(ceil(s));     % 'pseudo vdisc'
psvd = psvd(find(psvd<lenu));

% Create delayed unit with 1 sec. delay
flag = vdisc(find(vdisc>(lenu-10000)));
iflag = vdisc(find(vdisc<=(lenu-10000)));
dvd = [flag-lenu+10000 iflag+10000];     % 'delayed vdisc': last second of unit comes to first

% Create delayed unit with 10 sec. delay
flag2 = vdisc(find(vdisc>(lenu-100000)));
iflag2 = vdisc(find(vdisc<=(lenu-100000)));
dvd2 = [flag2-lenu+100000 iflag2+100000];     % 'delayed vdisc': last second of unit comes to first

% Wavelet transformation
dat = eeg(1:10:end);            % eeg
wavea = eeg_wavelet(dat);

dat = rand(size(dat));          % random
wavec = eeg_wavelet(dat);

% Sinc convolution & wavelet transformation
waveb = unit_wavelet(vdisc,lenu);     % unit
waved = unit_wavelet(psvd,lenu);     % random unit
wavee = unit_wavelet(dvd,lenu);     % delayed unit (1 sec. delay)
wavef = unit_wavelet(dvd2,lenu);     % delayed unit (10 sec. delay)



% --------------------------------------------------------------------------------------------
function wave = eeg_wavelet(dat)

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
mif = 0.5;          %minimal intresting frequency
mis = find(f>mif);
mis = mis(end);     %maximal intristing scale
mother = 'Morlet';

% Wavelet transformation
[wave,period,scale,coi] = b_wavelet_pow3(dat,dt,pad,dj,s0,j1,mother,param,mis);



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