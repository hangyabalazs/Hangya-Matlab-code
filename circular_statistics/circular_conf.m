function [L, U] = circular_conf(ang,degrad)
%CIRCULAR_CONF   Confidence interval for circular mean.
%   [L, U] = WCONF(ANG,DEGRAD) calculates (L,U) 95% confidence interval
%   from phase angles in ANG. The second input argument (DEGRAD) should contain
%   the information whether the sample is measured in degress ('deg') or
%   radians ('rad'). L and U are returned in the same measure.
%
%   See also CIRCULAR_MEAN3, MVL, MVLMN and CIRCULAR_SE.

% Input argument check
error(nargchk(2,2,nargin))

% Main
if isequal(degrad,'deg')    % convert to radian
    ang = ang / 180 * pi;
end
n = length(ang);    % sample size
alpha = 0.05;       % 95% conf. int.

[ftm, cm, R] = mvlmn(ang,'rad');    % first trignometric moment, circular mean and mean resultant length

C2 = 1 / n * sum(cos(2*ang));
S2 = 1 / n * sum(sin(2*ang));
m2 = sqrt(C2^2+S2^2);  % second central trigonometric moment
% m2 = (1 / n) * sum(cos(2*(ang-cm)));
delta = (1 - m2) / (2 * R^2);   % circular dispersion
sigma_square = delta / n;   % circular standard error
sigma = sqrt(sigma_square);
z = norminv(1-alpha/2,0,1);
L = cm - asin(z*sigma);     % confidence interval for mean direction
U = cm + asin(z*sigma);

if ~isreal(L)   % conf. int. covers the entire circle for very low mean resultant length
    L = NaN;
    U = NaN;
end

if isequal(degrad,'deg')    % convert to degrees
    L = L * 180 / pi;
    U = U * 180 / pi;
end