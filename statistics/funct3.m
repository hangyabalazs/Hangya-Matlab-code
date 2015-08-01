function fc = funct3(xc)
%FUNCT3   Error function for WATSONTWOFIT2.
%   FC = FUNCT3(XC) returns least square error in FC. Input argument XC
%   should be 1-by-4 array containing mu1, mu2, kappa1, kappa2 - mean and
%   concentration parameter of the two von Mises distributions.
%
%   FMINSEARCHBND calls FUNT3 by minimization.
%
%   See also FMINSEARCHBND and WATSONTWOFIT2.

% Get global variable 'THETA'
global THETA
theta = THETA;

% Define parameters
mu1 = xc(1);
mu2 = xc(2);
kappa1 = xc(3);
kappa2 = xc(4);

% Fix 'p'
p = 0.5;

% Calculate 'Rsquare'
n = length(theta);

A1_kappa1 = b_A1(kappa1);
A2_kappa1 = 1 - 2 / kappa1 * A1_kappa1;
A3_kappa1 = A1_kappa1 - 4 / kappa1 * A2_kappa1;
A1_kappa2 = b_A1(kappa2);
A2_kappa2 = 1 - 2 / kappa2 * A1_kappa2;
A3_kappa2 = A1_kappa2 - 4 / kappa2 * A2_kappa2;

C1 = sum(cos(theta)/n);
S1 = sum(sin(theta)/n);
C2 = sum(cos(2*theta)/n);
S2 = sum(sin(2*theta)/n);
C3 = sum(cos(3*theta)/n);
S3 = sum(sin(3*theta)/n);

deltaC1 = p * A1_kappa1 * cos(mu1) + (1 - p) * A1_kappa2 * cos(mu2) - C1;
deltaC2 = p * A2_kappa1 * cos(2*mu1) + (1 - p) * A2_kappa2 * cos(2*mu2) - C2;
deltaC3 = p * A3_kappa1 * cos(3*mu1) + (1 - p) * A3_kappa2 * cos(3*mu2) - C3;
deltaS1 = p * A1_kappa1 * sin(mu1) + (1 - p) * A1_kappa2 * sin(mu2) - S1;
deltaS2 = p * A2_kappa1 * sin(2*mu1) + (1 - p) * A2_kappa2 * sin(2*mu2) - S2;
deltaS3 = p * A3_kappa1 * sin(3*mu1) + (1 - p) * A3_kappa2 * sin(3*mu2) - S3;

Rsquare = deltaC1^2 + deltaC2^2 + deltaC3^2 + deltaS1^2 + deltaS2^2 + deltaS3^2;

% Output argument
fc = Rsquare;