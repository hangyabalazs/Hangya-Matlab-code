function [mu,kappa,Value,p,Rsquare] = watsontest(x)
%WATSONTEST     Watson's goodness of fit test for the von Mises distribution.
%   WATSONTEST(X) performs Watson's U^2-test on a vector of angular 
%   measurements given in radians (X).
%
%   WATSONTEST estimates the parameters of the von Mises distribution (mu
%   and kappa) and tests whether the given sample shows von Mises
%   distribution (goodness of fit test).
%
%   [MU,KAPPA,U2,P,R] = WATSONTEST(X) returns the estimated parameters (MU
%   and KAPPA), the test statistic (U2), limiting values for the p-value
%   (P) and the error of fit (R).
%
%   See also WATSONTWO and RAO.

%   Some parts of this code are modified from the circular statistics
%   package of R (www.r-rpoject.org).

%   Balazs Hangya
%   Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   New York 11724, USA
%   balazs.cshl@gmail.com
%
%   Institute of Experimental Medicine
%   Szigony street 43, Budapest 1083, Hungary
%   hangyab@koki.hu

% Input argument check
error(nargchk(1,1,nargin))

% Critical values
n = length(x);
u21 = [0, 0.5, 1, 1.5, 2, 4, 100]';
u22 = [0.052, 0.056, 0.066, 0.077, 0.084, 0.093, 0.096]';
u23 = [0.061, 0.066, 0.079, 0.092, 0.101, 0.113, 0.117]';
u24 = [0.081, 0.09, 0.11, 0.128, 0.142, 0.158, 0.164]';
u2_crits = [u21 u22 u23 u24];

% Parameter estimation
[ftm mu] = mvlmn(x);   % circular mean
kappa = A1inv(mean(cos(x-mu)));   % concentration parameter

% Watson's U^2-test
x = mod((x-mu),2*pi);
x = x(:);
z = pvm(x,0,kappa);   % von Mises cdf
z = sort(z);
zbar = mean(z);
i = 1:n;
sumterms = (z' - (2 * i - 1) / (2 * n)).^2;
Value = sum(sumterms) - n * (zbar - 0.5)^2 + 1 / (12 * n);   % U^2 stat.
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

if Value < u2_crits(row,2)   % p value
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
A1_kappa1 = A1(kappa);
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

Rsquare = deltaC1^2 + deltaC2^2 + deltaC3^2 + deltaS1^2 + deltaS2^2 + deltaS3^2;   % error

% -------------------------------------------------------------------------
function kappa = A1inv(x)

% Inverse function of the ratio of the first and zeroth order Bessel
% functions of the first kind.  This function is used to compute the
% maximum likelihood estimate of the concentration parameter of a
% von Mises distribution.
X = (0 <= x & x < 0.53);
YES = 2 * x + x^3 + (5 * x^5) / 6;
NO = ifelse(x<0.85,-0.4+1.39*x+0.43/(1-x),1/(x^3-4*x^2+3*x));
kappa = ifelse(X,YES,NO);

% -------------------------------------------------------------------------
function [ftm, ang, mvl] = mvlmn(ng,degrad)
%MVLMN   Circular mean and mean resultant length.
%   [F, M, R] = MVLMN(NG,DEGRAD) calculates the first trigonometric moment
%   (F), the circular mean (M) and the mean resultant length (R) of a given
%   circular sample (NG). The second input argument (DEGRAD) should contain
%   the information whether the sample is measured in degress ('deg') or
%   radians ('rad'). M is returned in the same measure.
%
%   See also MVL, CIRCULAR_MEAN3 and CIRCULAR_SE.

% Input argument check
error(nargchk(2,2,nargin))

% Main
if isequal(degrad,'deg')    % convert to radian
    ng = ng / 180 * pi;
end
ftm = sum(exp(1).^(i*ng)) / length(ng);    % first trigonometric moment
ang = angle(ftm);   % mean angle
mvl = abs(ftm);     % mean resultant length
if isequal(degrad,'deg')    % convert to degrees
    ang = ang * 180 / pi;
end

% -------------------------------------------------------------------------
function result = pvm(theta,mu,kappa,acc)
%PVM    Estimates the cummulative probability for a von Mises distribution.
%   P = PVM(THETA,MU,KAPPA,ACC) estimates the cummulative probability for a
%   von Mises distribution, where the input argumnets are:
%       theta: angular value in radians,
%       mu: mean direction of the von Mises distribution,
%       kappa: concentration parameter of the von Mises distribution,
%       acc: parameter related to the accuracy of the estimated 
%           cummulative probability.  See details below.  Default value
%           is 1e-020.
%
%   Cummulative probabilities are computed according to the expression for
%   the von Mises cdf given by Gumbel et al. (1953), which gives the cdf as
%   a function of an infinite sum.  The parameter acc specifies the 
%   accuracy with which this sum is approximated. Terms greater than acc
%   are included in the summation.
%
%   PVM returns the probability (P) that a von Mises random variable falls
%   between 0 and theta.
%
%   See also WATSON.

% Input argument check
error(nargchk(3,4,nargin))
if nargin == 3
    acc = 1e-020;
end

% Probability calculation
theta = mod(theta,2*pi);
mu = mod(mu,2*pi);
if mu == 0
    result = mu0(theta, kappa, acc);
else
    if theta <= mu
        upper = mod(theta-mu,2*pi);
        if upper == 0
            upper = 2 * pi;
        end
        lower = mod(-mu,2*pi);
        result = mu0(upper,kappa,acc) - mu0(lower,kappa,acc);
    else
        upper = theta - mu;
        lower = mod(mu,2*pi);
        result = mu0(upper,kappa,acc) + mu0(lower,kappa,acc);
    end
end

% -------------------------------------------------------------------------
function M = mu0(theta,kappa,acc)
flag = 1;
p = 1;
sm = 0;
while flag
    term = besseli(p,kappa) * sin(p*theta) / p;
    sm = sm + term;
    p = p + 1;
    if abs(term) < acc
        flag = 0;
    end
end
M = theta / (2 * pi) + sm / (pi * besseli(0,kappa));

% -------------------------------------------------------------------------
function a1 = A1(kappa)
%A1   Ratio of first and zeroth order Bessel functions.
%   A1(KAPPA) Evaluates the first and zeroth order Bessel functions of the
%   first kind at KAPPA non-negative real number, and returns the ratio.
%
%   See also BESSELI.

% Input argument check
error(nargchk(1,1,nargin))

% A1
a1 = besseli(1,kappa) ./ besseli(0,kappa); 

% -------------------------------------------------------------------------
function C = ifelse(x,yes,no)
%IFELSE    Conditional mixing.
%   IFELSE(X,YES,NO) mixes the elements of YES and NO in the following way:
%   it takes the element of YES if corresponding element of logical input X
%   is true, otherwise it takes the corresponding element of NO. It returns
%   an X-size matrix.
%
%   Note that YES and NO should be either X-size, or 1x1!
%
%   See also A1INV.

% Input argument check
error(nargchk(3,3,nargin)) 
if ~isequal(size(x),size(yes))
    if ~isequal(size(yes),[1 1])
        error('Input argumnets X and YES must be of equal size or YES should be 1-by-1.')
    else
        yes = repmat(yes,size(x));
    end
end
if ~isequal(size(x),size(no))
    if ~isequal(size(no),[1 1])
        error('Input argumnets X and NO must be of equal size or NO should be 1-by-1.')
    else
        no = repmat(no,size(x));
    end
end
    
% Mixing
A = x .* yes;
B = ~x .* no;
C = A + B;