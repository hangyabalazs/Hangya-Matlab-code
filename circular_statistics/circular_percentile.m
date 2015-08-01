function [pbl pab] = circular_percentile(ang,prc,degrad)
%CIRCULAR_PERCENTILE   Circular percentile.
%   [P1 P2] = CIRCULAR_PERCENTILE(ANG,PRC,DEGRAD) calculates circular
%   percentiles (PRC) below and above the mean phase from phase angles in
%   ANG. The third input argument (DEGRAD) should contain the information
%   whether the sample is measured in degress ('deg') or radians ('rad').
%   Lower (P1) and upper percentile (P2) values are returned in the same
%   measure.
%
%   See also CIRCULAR_MEAN3,, CIRCULAR_SE, MVL and MVLMN.

%   Balazs Hangya
%   Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   New York 11724, USA
%   balazs.cshl@gmail.com
%
%   Institute of Experimental Medicine
%   Szigony street 43, Budapest 1083, Hungary
%   hangyab@koki.hu

% Input argument check
error(nargchk(3,3,nargin))

% Conversion
if isequal(degrad,'deg')    % convert to radian
    ang = ang / 180 * pi;
end
ang = mod(ang,2*pi);    % values between 0 and 2*pi
n = length(ang);    % sample size

% First trig. moment
[ftm, cm, R] = mvlmn(ang,'rad');    % first trignometric moment, circular mean and mean resultant length
cm = mod(cm,2*pi);

% Upper percentile
angsabove = [sort(ang(ang>=cm),'ascend') sort(ang(ang<cm),'ascend')+2*pi];  % data points 'upwards' from mean
nmpab = n * prc / 100;   % number of data points within the upper percentile
nmpab1 = floor(nmpab);  % interpolate
nmpab2 = ceil(nmpab);
if isequal(nmpab1,nmpab2)   % integer number of data points, no need to interpolate
    mltp = 0;
else
    mltp = (nmpab - nmpab1) / (nmpab2 - nmpab1);
end
pab = angsabove(nmpab1) + mltp * (angsabove(nmpab2) - angsabove(nmpab1));  % upper percentile value

% Lower percentile
angsbelow = [sort(ang(ang<cm),'descend') sort(ang(ang>=cm),'descend')-2*pi];  % data points 'below' mean
pbl = angsbelow(nmpab1) + mltp * (angsbelow(nmpab2) - angsbelow(nmpab1));  % upper percentile value

% Conversion
pab = mod(pab,2*pi);
pab(pab>pi) = pab(pab>pi) - 2 * pi;   % output between -pi and pi by convention
pbl = mod(pbl,2*pi);
pbl(pbl>pi) = pbl(pbl>pi) - 2 * pi;   % output between -pi and pi by convention
if isequal(degrad,'deg')    % convert to degrees
    pab = pab * 180 / pi;
    pbl = pbl * 180 / pi;
end