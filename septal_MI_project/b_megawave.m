function [wavea,waveb] = b_megawave
%MEGAWAVE   Calculates wavelet transform for unit and eeg.
%   Unit should be discriminated first. The time series are preprocessed through sinc 
%   convolution before wavelet transformation.
%
%   [A,B] = MEGAWAVE returns two output arguments: eeg wavelet (A) and unit wavelet (B).
%
%   See also MEGAWAVERUN and MEGAWAVERUN2.

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

% Sinc convolution
fs = 10000;
dto = 1 / fs;
ts = zeros(1,length(unit));
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
tnew = (1:length(unit)*fsratio) * dtnew * fcut;
lentold = length(told);
zint = 0;
for i = 1:lentold
    zint = zint + sinc(tnew-told(i));
   if  mod(i,100) == 0
       disp('%')
   end
end

% Wavelet transformation
dat = eeg(1:10:end);
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
disp(' COMPUTING THE EEG WAVELET SPECTRUM ...')
[wavea,period,scale,coi] = b_wavelet_pow3(dat,dt,pad,dj,s0,j1,mother,param,mis);

dat = zint;
disp(' COMPUTING THE UNIT WAVELET SPECTRUM ...')
[waveb,period,scale,coi] = b_wavelet_pow3(dat,dt,pad,dj,s0,j1,mother,param,mis);