function kappa = A1inv(x)
%A1INV   Inverse function of the ratio of the first and zeroth order Bessel functions.
%   K = A1INV(X) calculates the inverse function of the ratio of the first 
%   and zeroth order Bessel functions of the first kind at X.  This
%   function is used to compute the maximum likelihood estimate of the
%   concentration parameter of a von Mises distribution. In the latter
%   case, the mean resultant length of the von Mises sample should be given
%   as input argument.
%
%   See also A1 and KAPPACOMPARE.

X = (0 <= x & x < 0.53);
YES = 2 * x + x^3 + (5 * x^5) / 6;
NO = b_ifelse(x<0.85,-0.4+1.39*x+0.43/(1-x),1/(x^3-4*x^2+3*x));
kappa = b_ifelse(X,YES,NO);