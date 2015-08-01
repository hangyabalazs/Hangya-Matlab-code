function sigma = circular_SE(ang,degrad)
%CIRCULAR_SE   Circular standard error.
%   S = CIRCULAR_SE(ANG, DEGRAD) calculates circular standard error from
%   phase angles in ANG. The second input argument (DEGRAD) should contain
%   the information whether the sample is measured in degress ('deg') or
%   radians ('rad'). S is returned in the same measure.
%
%   See also CIRCULAR_MEAN3, MVL and MVLMN.

% Input argument check
error(nargchk(2,2,nargin))

% Main
if isequal(degrad,'deg')    % convert to radian
    ang = ang / 180 * pi;
end
n = length(ang);    % sample size

[ftm, cm, R] = mvlmn(ang,'rad');    % first trignometric moment, circular mean and mean resultant length

C2 = 1 / n * sum(cos(2*ang));
S2 = 1 / n * sum(sin(2*ang));
m2 = sqrt(C2^2+S2^2);  % second central trigonometric moment
% m2 = (1 / n) * sum(cos(2*(ang-cm)));
delta = (1 - m2) / (2 * R^2);   % circular dispersion
sigma_square = delta / n;   % circular standard error
sigma = sqrt(sigma_square);
if isequal(degrad,'deg')    % convert to degrees
    sigma = sigma * 180 / pi;
end