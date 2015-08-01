function [hang_tc hmvl_tc hang_ct hmvl_ct cang_tc cmvl_tc cang_ct cmvl_ct]...
    = b_wangle(tI,cI,vdtI,vdcI,wave_cross)
%WANGLE   Hilbert and crosswavelet mean angle and mean vector length for Tisza data.
%   WANGLE(TI,CI,VDTI,VDCI,WAVE_CROSS) computes and returns Hilbert mean
%   angle and mean resultant length: angles are the values of 
%   Hilbert-transformated TI at locations of events in VDCI, or values of
%   Hilbert-transformated CI at the values of VDTI. It calculates 
%   crosswavelet angle: finds maximums in cross power of CI and TI, and
%   takes the angle values corresponding to the maximums at locations of
%   events (VDCI or VDTI). Input argument WAVE_CROSS should be the
%   crosswavelet of TI and CI, which can be calculated using WCROSSWAVELET.
%
%   See also WCROSSWAVELET and ZSHIFTRUN2.

[hang_tc hmvl_tc] = hangle(vdtI,cI);
[hang_ct hmvl_ct] = hangle(vdcI,tI);
[cang_tc cmvl_tc] = cangle(vdtI,wave_cross);
[cang_ct cmvl_ct] = cangle(vdcI,wave_cross);

% -------------------------------------------------------------------------
function [hang hmvl] = hangle(vdisc,eeg)

% Hilbert transformation of the eeg
ahee = angle(hilbert(eeg));
aheedeg = ahee * (180 / pi);

% Phase angles - Hilbert
bahee = ahee(vdisc);
n = length(bahee);
ftm = sum(exp(1).^(i*bahee)) / n;    % first trigonometric moment
hang = angle(ftm);   % mean angle
hmvl = abs(ftm);     % mean resultant length

% -------------------------------------------------------------------------
function [cang cmvl] = cangle(vdisc,wave_cross)

% Find band
pwind1 = 60;
pwind2 = 90;

% Phase angles - crosswavelet
wph = angle(wave_cross(pwind1:pwind2,:));
wabs = abs(wave_cross(pwind1:pwind2,:));
mwabs = max(wabs);
lenw = length(wabs);

wph0 = zeros(1,lenw);
for k = 1:lenw
    mloc = find(wabs(:,k)==mwabs(k));
    wph0(k) = wph(mloc,k);
end
wphh = wph0(vdisc);

n = length(wphh);
ftm = sum(exp(1).^(i*wphh)) / n;    % first trigonometric moment
cang = angle(ftm);   % mean angle
cmvl = abs(ftm);     % mean resultant length