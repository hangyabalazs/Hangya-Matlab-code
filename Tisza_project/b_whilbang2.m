function [hang_tc hmvl_tc hang_ct hmvl_ct] = b_whilbang2(tI,cI,vdtI,vdcI)
%WHILBANG2   Hilbert mean angle and mean vector length for Tisza data.
%   WHILBANG2(TI,CI,VDTI,VDCI) computes and returns Hilbert mean angle and
%   mean resultant length: angles are the values of filtered and Hilbert-
%   transformated TI at locations of events in VDCI, or values of filtered
%   and Hilbert-transformated CI at the values of VDTI.
%
%   WHILBANG2 filters the signal sequentially in week resolution.
%
%   See also WANGLE, WANGLE2 and WHILBANG.

hang_tc = [];
hmvl_tc = [];
hang_ct = [];
hmvl_ct = [];
for k = 1:52    % weekly
    fr1 = 52 / k;
    fr2 = 52 / max((k - 1),0.2857);     % 2 days if k == 0
    [hang_tc(end+1) hmvl_tc(end+1)] = hangle(vdtI,cI,fr1,fr2);
    [hang_ct(end+1) hmvl_ct(end+1)] = hangle(vdcI,tI,fr1,fr2);
end

% -------------------------------------------------------------------------
function [hang hmvl] = hangle(vdisc,eeg,fr1,fr2)

% Filtering
flt = fir1(512,[fr1/(365/2) fr2/(365/2)]);
feeg = filtfilt(flt,1,eeg);
% figure
% subplot(2,1,1)
% plot(eeg)
% subplot(2,1,2)
% plot(feeg)

% Hilbert transformation of the eeg
ahee = angle(hilbert(feeg));
aheedeg = ahee * (180 / pi);

% Phase angles - Hilbert
bahee = ahee(vdisc);
n = length(bahee);
ftm = sum(exp(1).^(i*bahee)) / n;    % first trigonometric moment
hang = angle(ftm);   % mean angle
hmvl = abs(ftm);     % mean resultant length