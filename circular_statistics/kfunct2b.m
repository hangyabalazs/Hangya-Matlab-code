function fc = kfunct2b(xc)
%KFUNCT2B   Error function for WATSONKFIT.
%   FC = KFUNCT2B(XC) returns least square error in FC. Input argument XC
%   should be 1-by-3K-1 array containing circular means and concentration
%   parameters of the von Mises distributions and mixing ratios.
%
%   FMINSEARCHBND calls KFUNT2B by minimization.
%
%   See also FMINSEARCHBND and WATSONKFIT.

% Get global variable 'THETA'
global THETA
theta = THETA;

% Define parameters
k = (length(xc) + 1) / 3;
mu = xc(1:k);
kappa = xc(k+1:2*k);
p = xc(2*k+1:end);
p(end+1) = 1 - sum(p);

% Calculate 'Rsquare'
n = length(theta);

A1_kappa = b_A1(kappa);
A2_kappa = 1 - 2 ./ kappa .* A1_kappa;
A3_kappa = A1_kappa - 4 ./ kappa .* A2_kappa;

C1 = sum(cos(theta)/n);
S1 = sum(sin(theta)/n);
C2 = sum(cos(2*theta)/n);
S2 = sum(sin(2*theta)/n);
C3 = sum(cos(3*theta)/n);
S3 = sum(sin(3*theta)/n);

deltaC1 = sum(p.*A1_kappa.*cos(mu)) - C1;
deltaC2 = sum(p.*A2_kappa.*cos(2*mu)) - C2;
deltaC3 = sum(p.*A3_kappa.*cos(3*mu)) - C3;
deltaS1 = sum(p.*A1_kappa.*sin(mu)) - S1;
deltaS2 = sum(p.*A2_kappa.*sin(2*mu)) - S2;
deltaS3 = sum(p.*A3_kappa.*sin(3*mu)) - S3;

Rsquare = deltaC1^2 + deltaC2^2 + deltaC3^2 + deltaS1^2 + deltaS2^2 + deltaS3^2;

% Output argument
fc = Rsquare;