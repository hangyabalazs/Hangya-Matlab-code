function [mu,kappa,Value,p,Rsquare] = b_watson2(x,pl)
%WATSON2     Watson's goodness of fit test for the von Mises distribution.
%   WATSON2(X) performs Watson's U^2-test on X vector of angular 
%   measurements given in radians.
%
%   WATSON2 estimates the parameters of the von Mises distribution (mu and
%   kappa) and tests whether the given sample shows von Mises distribution
%   (goodness of fit test).
%
%   [MU,KAPPA,U2,P,R] = WATSON2(X) returns the estimated parameters (MU and
%   KAPPA), the test statistic (U2), limiting values for the p-value (P)
%   and the error of fit (R).
%
%   WATSON2(X,PL) plots the result of the fit, if PL = 1.
%
%   See also WATSONTWO and RAO.

% Input argument check
error(nargchk(2,2,nargin))

% Critical values
n = length(x);
u21 = [0, 0.5, 1, 1.5, 2, 4, 100]';
u22 = [0.052, 0.056, 0.066, 0.077, 0.084, 0.093, 0.096]';
u23 = [0.061, 0.066, 0.079, 0.092, 0.101, 0.113, 0.117]';
u24 = [0.081, 0.09, 0.11, 0.128, 0.142, 0.158, 0.164]';
u2_crits = [u21 u22 u23 u24];

% Parameter estimation
mu = b_circular_mean3(x);
kappa = A1inv(mean(cos(x-mu)));

% Watson's U^2-test
theta = mod(x,2*pi);
x = mod((x-mu),2*pi);
x = x(:);
z = b_pvm(x,0,kappa);
z = sort(z);
zbar = mean(z);
i = [1:n];
sumterms = (z' - (2 * i - 1) / (2 * n)).^2;
Value = sum(sumterms) - n * (zbar - 0.5)^2 + 1 / (12 * n);
if kappa < 0.25
    row = 1;
elseif kappa < 0.75
    row = 2;
elseif kappa < 1.25
    row = 3;
elseif kappa < 1.75
    row = 4;
elseif kappa < 3
    row = 5;
elseif kappa < 5
    row = 6;
else
    row = 7;
end

if Value < u2_crits(row,2)
    p(1) = 0.10;
    p(2) = 1;
elseif  (Value >= u2_crits(row,2)) && (Value < u2_crits(row,3))
    p(1) = 0.05;
    p(2) = 0.10;
elseif (Value >= u2_crits(row,3)) && (Value < u2_crits(row,4))
    p(1) = 0.01;
    p(2) = 0.05;
else
    p(1) = 0;
    p(2) = 0.01;
end

% Error of fit (Rsquare)
A1_kappa1 = b_A1(kappa);
A2_kappa1 = 1 - 2 / kappa * A1_kappa1;
A3_kappa1 = A1_kappa1 - 4 / kappa * A2_kappa1;

C1 = sum(cos(x)/n);
S1 = sum(sin(x)/n);
C2 = sum(cos(2*x)/n);
S2 = sum(sin(2*x)/n);
C3 = sum(cos(3*x)/n);
S3 = sum(sin(3*x)/n);

deltaC1 = A1_kappa1 * cos(mu) - C1;
deltaC2 = A2_kappa1 * cos(2*mu) - C2;
deltaC3 = A3_kappa1 * cos(3*mu) - C3;
deltaS1 = A1_kappa1 * sin(mu) - S1;
deltaS2 = A2_kappa1 * sin(2*mu) - S2;
deltaS3 = A3_kappa1 * sin(3*mu) - S3;

Rsquare = deltaC1^2 + deltaC2^2 + deltaC3^2 + deltaS1^2 + deltaS2^2 + deltaS3^2;

if pl
%     figure
    xx = linspace(0,2*pi,length(x));
    dvm = 1 / (2 * pi * besseli(0,kappa,1)) * (exp(cos(xx-mu)-1)).^kappa;
    hold on
    bno = fix(exp(0.626+0.4*log(length(x)-1)));   % number of bins
    [thb,nb] = hist(theta,bno);
    thb_norm = thb / mean(thb) * mean(dvm);
    bar([nb nb+2*pi],[thb_norm thb_norm]);
    plot([xx xx+2*pi],[dvm dvm],'r','LineWidth',2)
end

% -------------------------------------------------------------------------
function kappa = A1inv(x)

% Inverse function of the ratio of the first and zeroth order Bessel
% functions of the first kind.  This function is used to compute the
% maximum likelihood estimate of the concentration parameter of a
% von Mises distribution.
X = (0 <= x & x < 0.53);
YES = 2 * x + x.^3 + (5 * x.^5) / 6;
NO = b_ifelse(x<0.85,-0.4+1.39*x+0.43./(1-x),1./(x.^3-4*x.^2+3*x));
kappa = b_ifelse(X,YES,NO);