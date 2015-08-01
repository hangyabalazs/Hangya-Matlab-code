function [Y p] = circmeancomp(ng1,ng2,degrad)
%CIRCMEANCOMP    Test for a common mean direction of two circular samples.
%   [Y,P] = CIRCMEANCOMP(X,Y,DEGRAD) performs test for a common circular
%   mean of two samples (X and Y) measured in radians (DEGRAD = 'rad') or
%   degrees (DEGRAD = 'deg'). Test statistic (Y) and p-value (P) is
%   returned. The test is appropriate if each sample size is larger than
%   25, and the ratio of sample circular dispersions does not exceed 4.
%
%   Reference:
%   Fisher, N. (1993) Statistical Analysis of Circular Data. Cambridge:
%   Cambridge Univerity Press, pp. 115-117 ('Method P')
%   
%   See also WATSONTWO and RAO.

% Input argument check
error(nargchk(3,3,nargin))
if isequal(degrad,'deg')    % convert to radian
    ng1 = ng1 / 180 * pi;
    ng2 = ng2 / 180 * pi;
end

% C, S, R
n1 = length(ng1);
n2 = length(ng2);
N = n1 + n2;
[ftm1, mu1, R1] = mvlmn(ng1,'rad');    % first trignometric moment, circular mean and mean resultant length
[ftm2, mu2, R2] = mvlmn(ng2,'rad');
C = n1 * cos(mu1) + n2 * cos(mu2);
S = n1 * sin(mu1) + n2 * sin(mu2);
R = sqrt(C^2 + S^2);

% Circular dispersion 
delta1 = circular_dispersion(ng1,'rad');   % circular dispersion
delta2 = circular_dispersion(ng2,'rad');
delta0 = (n1 * delta1 + n2 * delta2) / N;
dr = max(delta1,delta2) / min(delta1,delta2);
if dr > 4
    warning('Ratio of circular dispersions exceed 4.')
end

% Test statistic
Y = 2 * (N - R) / delta0;

% Test
alpha = 0.05;
c = chi2inv(1-alpha,1);
p = 1 - chi2cdf(Y,1);