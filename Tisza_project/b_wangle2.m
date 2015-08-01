function [hang_tc hmvl_tc hang_ct hmvl_ct cang_tc cmvl_tc cang_ct cmvl_ct angvec_tc angvec_ct]...
    = b_wangle2(tI,cI,vdtI,vdcI,wave_cross)
%WANGLE2   Hilbert and crosswavelet mean angle and mean vector length for Tisza data.
%   WANGLE2(TI,CI,VDTI,VDCI,WAVE_CROSS) computes and returns Hilbert mean
%   angle and mean resultant length: angles are the values of 
%   Hilbert-transformated TI at locations of events in VDCI, or values of
%   Hilbert-transformated CI at the values of VDTI. It calculates 
%   crosswavelet angle: finds maximums in cross power of CI and TI, and
%   takes the angle values corresponding to the maximums at locations of
%   events (VDCI or VDTI). WANGLE2 returns crosswavelet angle vector, 
%   mean vector length and mean angle at the values of VDTI and VDCI. 
%   Input argument WAVE_CROSS should be the crosswavelet of TI and CI, 
%   which can be calculated using WCROSSWAVELET.
%
%   WANGLE2, in contrast to WANGLE calculates Hilbert angle and mean
%   resultant length from "different" angle values, i.e. if the difference
%   between consecutive values does not reach a threshold of 0.05, it
%   replaces them with their circular mean.
%
%   See also WANGLE, WCROSSWAVELET and ZSHIFTRUN2.

[hang_tc hmvl_tc] = hangle(vdtI,cI);
[hang_ct hmvl_ct] = hangle(vdcI,tI);
[cang_tc cmvl_tc angvec_tc] = cangle(vdtI,wave_cross);
[cang_ct cmvl_ct angvec_ct] = cangle(vdcI,wave_cross);

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
function [cang cmvl wphh] = cangle(vdisc,wave_cross)

% Find band
pwind1 = 60;
pwind2 = 100;

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

dwphh = diff(wphh);    % get phase "different" values
fdw = [0 abs(dwphh)>0.05];
% dfdw = diff(fdw);
% fdf1 = find([-1 dfdw]==-1);
% fdf2 = find(dfdw==1);
% dfdw = [-1 dfdw];
wnew = [];
for i = 1:length(fdw)
    switch fdw(i)
        case 0
            if exist('act')
                act(end+1) = wphh(i);
            else
                act = wphh(i);
            end
        case 1
            if i>1 && fdw(i-1)==0
                wnew(end+1) = cmean(act);
                act = [];
            elseif i>1 && fdw(i-1)==1
                wnew(end+1) = wphh(i-1);
            end
    end
end
wphh = wnew;

n = length(wphh);
ftm = sum(exp(1).^(j*wphh)) / n;    % first trigonometric moment
cang = angle(ftm);   % mean angle
cmvl = abs(ftm);     % mean resultant length

% -------------------------------------------------------------------------
function M = cmean(A)

ftm = sum(exp(1).^(j*A)) / length(A);    % first trigonometric moment
M = angle(ftm);   % mean angle