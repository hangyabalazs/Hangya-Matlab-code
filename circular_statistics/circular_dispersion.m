function delta = circular_dispersion(ang,degrad)
%CIRCULAR_DISPERSION   Circular dispersion.
%   DELTA = CIRCULAR_DISPERSION(NG,DEGRAD) calculates circular dispersion 
%   (DELTA) of a given circular sample (NG). The second input argument
%   (DEGRAD) should contain the information whether the sample is measured
%   in degress ('deg') or radians ('rad').
%
%   See also MVL, MVLMN, CIRCULAR_MEAN, CIRCULAR_CONF and CIRCULAR_SE.

% Input argument check
error(nargchk(2,2,nargin))
if isequal(degrad,'deg')    % convert to radian
    ang = ang / 180 * pi;
end

% Main
n = length(ang);    % sample size
[ftm, cm, R] = mvlmn(ang,'rad');    % first trignometric moment, circular mean and mean resultant length
C2 = 1 / n * sum(cos(2*ang));
S2 = 1 / n * sum(sin(2*ang));
m2 = sqrt(C2^2+S2^2);  % second central trigonometric moment
% m2 = (1 / n) * sum(cos(2*(ang-cm)));
delta = (1 - m2) / (2 * R^2);   % circular dispersion