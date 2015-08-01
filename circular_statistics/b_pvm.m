function result = b_pvm(theta,mu,kappa,acc)
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