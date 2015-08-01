function [wavea,waveb,wavec,waved] = b_megawave_entrctrl
%MEGAWAVE_ENTRCTRL   Calculates wavelet transform for unit, eeg, random sample and random generated unit.
%   Unit should be discriminated first. The time series are preprocessed through sinc 
%   convolution before wavelet transformation.
%
%   [A,B,C,D] = MEGAWAVE_ENTRCTRL returns 4 output arguments: eeg wavelet (A), unit wavelet (B),
%   random sample wavelet (C) and random generated unit wavelet (D). Random generated unit is a
%   randomized Poisson process.
%
%   See also MEGAWAVE_ENTRCTRL2 and MEGAWAVERUN_ENTRCTRL.

% Input arguments check
error(nargchk(0,0,nargin));

% Input variables
global IN
data = IN{1};
eeg = IN{2};
fname = IN{3};
pathname = IN{4};
datinx1 = IN{5};
datinx2 = IN{6};
time = IN{7};
unit = IN{8};
dt = IN{9};
meret = IN{10};
mintafr = IN{11};
xlimit = IN{12};

% Discrimination variables
global DISC
output = DISC{2};
vdisc = DISC{3};
kuszob = DISC{4};
instfrek = DISC{5};
isi = DISC{6};

% Create random unit
lenu = length(unit);
r = random('exp',1000,1,1500);
s(1) = r(1);
for i = 2:length(r)
    s(i) = s(i-1) + r(i);
end
r2 = cumsum(s);
psvd = ceil(s);
% z = [];
% z(psvd) = 1;
% z = z(1:lenu);
psvd = psvd(find(psvd<lenu));

% Wavelet transformation
dat = eeg(1:10:end);            % eeg
wavea = eeg_wavelet(dat);

dat = rand(size(dat));          % random
wavec = eeg_wavelet(dat);

% Sinc convolution & wavelet transformation
waveb = unit_wavelet(vdisc,lenu);     % unit
waved = unit_wavelet(psvd,lenu);     % unit



% --------------------------------------------------------------------------------------------
function wave = eeg_wavelet(dat)

% Prepare for wavelet transformation
variance = std(dat) ^ 2;
dat = (dat - mean(dat)) / sqrt(variance) ;
n = length(dat);
dt = 1 / 1000;
pad = 1;
dj = 0.08;    
j1 = ceil((1/dj)*log2(n/2));
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
disp(' COMPUTING THE EEG WAVELET SPECTRUM ...')     % eeg
[wave,period,scale,coi] = b_wavelet_new2(dat,dt,pad,dj,s0,j1,mother,param,mis);



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
   if  mod(i,100) == 0
       disp('%')
   end
end

% Prepare for wavelet transformation
variance = std(zint) ^ 2;
zint = (zint - mean(zint)) / sqrt(variance) ;
n = length(zint);
dt = 1 / 1000;
pad = 1;
dj = 0.08;    
j1 = ceil((1/dj)*log2(n/2));
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
disp(' COMPUTING THE EEG WAVELET SPECTRUM ...')
[wave,period,scale,coi] = b_wavelet_new2(zint,dt,pad,dj,s0,j1,mother,param,mis);