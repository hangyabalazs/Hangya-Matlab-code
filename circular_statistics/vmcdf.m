function F = vmcdf(theta,mu,kappa,degrad,tol)
%VMCDF   Distribution function for the von Mises distribution.
%   F = VMCDF(THETA,MU,KAPPA,DEGRAD,TOL) calculates approximate
%   distribution function (cdf) for the von Mises distribution
%   VM(MU,KAPPA) at values in THETA. The input argument (DEGRAD) should
%   contain the information whether the sample is measured in degress
%   ('deg') or radians ('rad'). TOL is the tolerance level of the
%   approximation (default: e-20). Values of the distribution function are
%   returned in F. See the reference below for mathematical detals.
%
%   Reference: Algorithm AS 86: The Von Mises Distribution Function
%   K. V. Mardia and P. J. Zemroch; Journal of the Royal Statistical 
%   Society. Series C (Applied Statistics), Vol. 24, No. 2 (1975), pp. 
%   268-272
%
%   See also VMCOMPONENTS.

% Input argument check
error(nargchk(4,5,nargin))
if nargin < 5   % tolerance for the approximation
    tol = 10 ^(-20);
end
if isequal(degrad,'deg')    % convert to radian
    theta = theta / 180 * pi;
end

% Approximation
theta = [theta-mu -mu];
n = length(theta);

Q = zeros(1,n);
ad = ones(1,n);
k = 1;
while any(abs(ad)>tol)
    ad = besseli(k,kappa) * sin(k*theta) / k;
    
    Q = Q + ad;
    k = k + 1;
end

% Distribution function
F0 = 1 / (2 * pi * besseli(0,kappa)) * (besseli(0,kappa) * theta + 2 * Q);
F = F0(1:end-1) - F0(end);