function [kappa_ML kappa_hat] = kappaest(ng,degrad)
%KAPPAEST   Concentration parameter estimation for the von Mises distribution.
%   [KML, KHAT] = KAPPAEST(NG,DEGRAD) estimates the concentration parameter
%   of the von Mises sample NG. The second input argument (DEGRAD) should
%   contain the information whether the sample is measured in degress
%   ('deg') or radians ('rad'). Maximum likelihood estimate (KML) and
%   corrected maximum likelihood estimate (KHAT) are returned.
%
%   See also WATSON and KAPPACOMPARE.

% Input argument check
error(nargchk(2,2,nargin))

% Mean resultant length
n = length(ng);
[ftm, ang, mvl] = mvlmn(ng,degrad);

% Maximum likelihood estimate
if mvl < 0.53
    kappa_ML = 2 * mvl + mvl^3 + 5 * mvl^5 / 6;
elseif mvl < 0.85
    kappa_ML = -0.4 + 1.39 * mvl + 0.43 / (1 - mvl);
else
    kappa_ML = 1 / (mvl^3 - 4 * mvl^2 + 3 * mvl);
end
kappa_ML2 = A1inv(mvl);     % alternative way of calculation
if ~isequal(kappa_ML,kappa_ML2)
    error('Technical error 28')
end

% Corrected maximum likelihood estimate
if n <= 15
    if kappa_ML < 2
        kappa_hat = max(kappa_ML-2/(n*kappa_ML),0);
    else
        kappa_hat = (n-1)^3 * kappa_ML / (n^3 + n);
    end
else
    kappa_hat = kappa_ML;
end