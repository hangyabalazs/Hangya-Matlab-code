function [cang_tc cmvl_tc cang_ct cmvl_ct angvec_tc angvec_ct] = b_wcrossang(tI,cI,vdtI,vdcI,wave_cross)
%WCROSSANG   Crosswavelet mean angle and mean vector length for Tisza data.
%   WCROSSANG(TI,CI,VDTI,VDCI,WAVE_CROSS) computes and returns crosswavelet
%   mean angle and mean resultant length: angles are the phase values of 
%   wavelet-transformated TI at locations of events in VDCI, or values of
%   wavelet-transformated CI at the values of VDTI. It calculates 
%   crosswavelet angle as follows: finds maximums in cross power of CI and TI,
%   and takes the angle values corresponding to the maximums at locations of
%   events (VDCI or VDTI), in case corresponding maximum exceed a constant 
%   threshold. WCROSSANG returns crosswavelet angle vector, mean vector length
%   and mean angle at the values of VDTI and VDCI. Input argument WAVE_CROSS 
%   should be the crosswavelet of TI and CI, which can be calculated using 
%   WCROSSWAVELET.
%
%   See also WANGLE, WCROSSWAVELET and ZSHIFTRUN2.

[cang_tc cmvl_tc angvec_tc] = cangle(vdtI,wave_cross);
[cang_ct cmvl_ct angvec_ct] = cangle(vdcI,wave_cross);



% -------------------------------------------------------------------------
function [cang cmvl wphh] = cangle(vdisc,wave_cross)

% Find band
pwind1 = 20;
pwind2 = 30;

% Phase angles - crosswavelet
wph = angle(wave_cross(pwind1:pwind2,:));
wabs = abs(wave_cross(pwind1:pwind2,:)).^2;
mwabs = max(wabs);
lenw = length(wabs);

wph0 = zeros(1,lenw);
for k = 1:lenw
    mloc = find(wabs(:,k)==mwabs(k));
    wph0(k) = wph(mloc,k);
end
wphh = wph0(vdisc);
wahh = mwabs(vdisc);

% dwphh = diff(wphh);    % get "different" phase values
% fdw = [0 abs(dwphh)>0.05];
% wnew = [];
% for i = 1:length(fdw)
%     switch fdw(i)
%         case 0
%             if exist('act')
%                 act(end+1) = wphh(i);
%             else
%                 act = wphh(i);
%             end
%         case 1
%             if i>1 && fdw(i-1)==0
%                 wnew(end+1) = cmean(act);
%                 act = [];
%             elseif i>1 && fdw(i-1)==1
%                 wnew(end+1) = wphh(i-1);
%             end
%     end
% end
% wphh = wnew;

threshold = 10;     % drop subthreshold values
wphh = wphh(find(wahh>threshold));

n = length(wphh);
ftm = sum(exp(1).^(j*wphh)) / n;    % first trigonometric moment
cang = angle(ftm);   % mean angle
cmvl = abs(ftm);     % mean resultant length

% -------------------------------------------------------------------------
function M = cmean(A)

ftm = sum(exp(1).^(j*A)) / length(A);    % first trigonometric moment
M = angle(ftm);   % mean angle