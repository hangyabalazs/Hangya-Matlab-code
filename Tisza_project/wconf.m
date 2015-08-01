function [L, U] = wconf(ang,cm,R)
%WCONF   Confidence interval for circular mean.
%   [L, U] = WCONF(ANG,CM,R) calculates (L,U) 95% confidence interval for
%   circular mean (CM) from phase angles in ANG. R is the mean resultant 
%   length of ANG. Input arguments can be calculated via WCROSSANG. Note,
%   that L and U are given in radians!
%
%   See also WCROSSWAVELET, WCROSSWAVELET2 and WCROSSANG.

n = length(ang);    % sample size
alpha = 0.05;

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