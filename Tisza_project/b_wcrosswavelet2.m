function [wave_cross, f] = b_wcrosswavelet2(tI,cI)
%WCROSSWAVELET2   Crosswavelet for Tisza data.
%   WCROSSWAVELET2(TI,CI) calculates crosswavelet of TI and CI.
%   WCROSSWAVELET2 filters data instead of discrimination and convolution.
%
%   See also WAVELET_NEW3 and WCROSSWAVELET.

% Iput argument check
error(nargchk(2,2,nargin))

% Main
[wave_tI,f] = wavelet(tI,365);
[wave_cI,f] = wavelet(cI,365);
wave_cross = wave_tI .* conj(wave_cI);

% -------------------------------------------------------------------------
function [wave,f] = wavelet(eeg,sfr)

% Filtering
flt = fir1(512,[0.5/(365/2) 52/(365/2)]);        % 1 week - 2 year period time
feeg = filtfilt(flt,1,eeg);
feeg=eeg;

% Standardization
variance = std(feeg) ^ 2;
feeg = (feeg - mean(feeg)) / sqrt(variance) ;

% Prepare for wavelet transformation
n = length(feeg);
dt = 1 / sfr;
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
[wave,period,scale,coi] = b_wavelet_new3(feeg,dt,pad,dj,s0,j1,mother,param,mis);