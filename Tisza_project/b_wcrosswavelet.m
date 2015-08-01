function [wave_cross, f] = b_wcrosswavelet(vdtI,lentI,vdcI,lencI)
%WCROSSWAVELET   Crosswavelet for Tisza data.
%   WCROSSWAVELET(VDTI,LENTI,VDCI,LENCI) needs for input argument:
%       VDTI: discriminated TI (output of WDISC)
%       LENTI: length(TI)
%       VDCI: discriminated CI (output of WDISC)
%       LENCI: length(CI) (only for symmetry, equal to LENTI)
%
%   See also WAVELET_NEW3 and WCROSSWAVELET2.

% Iput argument check
error(nargchk(4,4,nargin))

% Main
[wave_tI,f] = wavelet(vdtI,lentI,365);
[wave_cI,f] = wavelet(vdcI,lencI,365);
wave_cross = wave_tI .* conj(wave_cI);

% -------------------------------------------------------------------------
function [wave,f] = wavelet(vdisc,lenu,sfr)

% Sinc convolution
fs = sfr;     % unit
dto = 1 / fs;
fcut = 100; 
fsnew = sfr;    % no downsampling
dtnew = 1 / fsnew;
fsold = fs;
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
dt = 1 / fsnew;
pad = 1;
dj = 0.1;
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
mis = -1;
% mif = 0.5;          %minimal intresting frequency
% mis = find(f>mif);
% mis = mis(end);     %maximal intristing scale
mother = 'Morlet';

% Wavelet transformation
[wave,period,scale,coi] = b_wavelet_new3(zint,dt,pad,dj,s0,j1,mother,param,mis);