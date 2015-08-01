function [hang_tc hmvl_tc hang_ct hmvl_ct] = b_whilbang(tI,cI,vdtI,vdcI)
%WHILBANG   Hilbert mean angle and mean vector length for Tisza data.
%   WHILBANG(TI,CI,VDTI,VDCI) computes and returns Hilbert mean angle and
%   mean resultant length: angles are the values of filtered and Hilbert-
%   transformated TI at locations of events in VDCI, or values of filtered
%   and Hilbert-transformated CI at the values of VDTI.
%
%   See also WANGLE and WANGLE2.

[hang_tc hmvl_tc] = hangle(vdtI,cI);
[hang_ct hmvl_ct] = hangle(vdcI,tI);

% -------------------------------------------------------------------------
function [hang hmvl] = hangle(vdisc,eeg)

% Filtering
flt = fir1(512,[1/(365/2) 2/(365/2)]);
feeg = filtfilt(flt,1,eeg);
figure
subplot(2,1,1)
plot(eeg)
subplot(2,1,2)
plot(feeg)

% Hilbert transformation of the eeg
ahee = angle(hilbert(feeg));
aheedeg = ahee * (180 / pi);

% Phase angles - Hilbert
bahee = ahee(vdisc);
n = length(bahee);
ftm = sum(exp(1).^(i*bahee)) / n;    % first trigonometric moment
hang = angle(ftm);   % mean angle
hmvl = abs(ftm);     % mean resultant length