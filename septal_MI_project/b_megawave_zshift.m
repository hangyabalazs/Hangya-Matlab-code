function [wavea,wavez] = b_megawave_zshift(eeg,vdisc,filenam)
%MEGAWAVE_ZSHIFT   Calculates wavelet transform for delayed unit.
%   EEG and 'vdisc' are needed as input arguments. The time series are 
%   preprocessed through sinc convolution before wavelet transformation.
%
%   [A,B] = MEGAWAVE_ENTRCTRL(EEG,VDISC) returns 2 output arguments: eeg 
%   wavelet (A) and shifted unit wavelet (B), where 'delay' is set to  
%   z-shift (see ZSHIFT for details).
%
%   See also MEGAWAVE_ENTRCTRL, ZSHIFT and MEGAWAVERUN_ZSHIFT.

% Input arguments check
error(nargchk(3,3,nargin));

% Create shifted unit
global DATAPATH
fn = [DATAPATH 'Zshift\Zshift3\text\zshift_theta3.xls'];
mtx = xlsread(fn);
[temp ln1] = xlsread(fn,'A2:A628');
inx = find(strcmp(filenam,ln1));
if length(inx) > 1
    zmx = mtx(inx,4);
    fnd = find(zmx==max(zmx));
    inx = inx(fnd);
end
dly = mtx(inx,3);
lenu = length(eeg);
dvd = vdisc(find(vdisc+dly>0&vdisc+dly<lenu)) + dly;  % shift unit

% Wavelet transformation
dat = eeg(1:10:end);            % eeg
wavea = eeg_wavelet(dat);

% Sinc convolution & wavelet transformation
wavez = unit_wavelet(dvd,lenu);     % zshift-delayed unit



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