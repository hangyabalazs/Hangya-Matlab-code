function x = vmpdf(theta,mu,kappa)
%VMPDF   Probability density function for the von Mises distribution.
%   X = VMPDF(THETA,MU,KAPPA) calculates the value of the von Mises pdf
%   with parameters MU and KAPPA at THETA.
%
%   See also RVM.

% Input argument check
error(nargchk(3,3,nargin))

% von Mises pdf
c = 1 / (2 * pi * besseli(0,kappa));
x = c * exp(kappa*cos(theta-mu));