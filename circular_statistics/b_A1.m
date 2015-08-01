function a1 = b_A1(kappa)
%A1   Ratio of first and zeroth order Bessel functions.
%   A1(KAPPA) Evaluates the first and zeroth order Bessel functions of the
%   first kind at KAPPA non-negative real number, and returns the ratio.
%
%   See also BESSELI.

% Input argument check
error(nargchk(1,1,nargin))

% A1
a1 = besseli(1,kappa) ./ besseli(0,kappa); 