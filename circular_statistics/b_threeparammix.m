function param = b_threeparammix(theta)
%THREEPARAMMIX   '3 parameter mix' fit of two von Mises mixture.
%   PARAM = THREEPARAMMIX(THETA) fits a mixture of two von Mises
%   distributions on the circular data set THETA. Note that THETA should
%   contain radian values! PARAM is a 1-by-5 array containing mu1, mu2,
%   kappa1, kappa2 - mean and concentration parameter of the two von Mises
%   distributions - and p mixing ratio (0 < p < 1).
%
%   See WATSON, WATSONTWO, WATSONTWOFIT, RVM and RMIXEDVM.

% Input argument check
error(nargchk(1,1,nargin))

% Estimate parameters
[mu0, kappa0, p0] = tpmx(theta);
mu1 = mu0;
mu2 = mu0 + pi;
kappa1 = kappa0;
kappa2 = kappa0;
p = p0;

% Output argument
param = [mu1 mu2 kappa1 kappa2 p];

% -------------------------------------------------------------------------
function [mu0, kappa0, p0] = tpmx(theta)

% Calculate accessoric variables
n = length(theta);
psi = mod(theta,pi);
Cpsi = sum(cos(2*psi));
Spsi = sum(sin(2*psi));
Rpsi = sqrt(Cpsi^2 + Spsi^2);
mupsi = b_circular_mean3(mod(theta,2*pi)-pi);
C1 = sum(cos(theta)/n);
S1 = sum(sin(theta)/n);

% Estimate 'mu'
mu0 = mupsi / 2;

% Method-of-moments estimate of 'kappa'
A1theta = b_A1(theta);  % linear interpolation for estimating inverse A2 function
A2theta = 1 - 2 ./ (theta + eps) .* A1theta;
fcn = abs(A2theta-Rpsi);
ind2 = find(fcn==min(fcn));
ind1 = ind2 - 1;
a = (ind2 - ind1) / (1 + fcn(ind2) / fcn(ind1));
x = ind1 + a;
kappa00 = A2theta(round(x));
                                                                                                     
% Estimate 'p0'
p0 = ((C1 * cos(mu0) + S1 * sin(mu0)) / b_A1(kappa00) + 1) / 2;
p0 = min(max(0,p0),1);      % 0 < p0 < 1

% Adjust 'kappa'
if kappa00 < 2  % 'kappa00 may be biased
    kappa0 = max(kappa00-2*(n*kappa00)^(-1),0);
else
    kappa0 = ((n - 1)^3 * kappa00) / (n^3 + n);
end
kappa0 = min(max(0,kappa0),1);      % 0 < kappa0 < 1