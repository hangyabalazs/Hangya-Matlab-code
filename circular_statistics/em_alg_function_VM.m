function [mu,kappa,Pa,iter,Q_tot,e_tot,P]=em_alg_function_VM(x,mu,kappa,Pa,bl,bu,e_min,maxiter)
%EM_ALG_FUNCTION_VM   Expectation Maximization for mixtures of von Mises distributions.
%   Input arguments:
%   x:  circular sample
%   mu: circular means for initial estimation
%   kappa: concenrtation parameters for initial estimation
%   Pa: initial mixing probabilities
%   bl: lower bound for the parameters
%   bu: upper bound for the parameters
%   e_min: tolerance (see below)
%   maxiter: maximum number of iterations
%
%   Outpu parameters:
%   mu: circular means
%   kappa: concentration parameters
%   Pa: mixing probabilities
%   iter: the number of iterations required for the convergence of the
%       EM algorithm.
%   Q_tot:  vector containing the likelihood value at each iteration.
%   e_tot:  vector containing the error value at each itertion.
%   P: conditional probabilities of the component from which the 
%   observations originates

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION
%   [m,s,Pa,iter,Q_tot,e_tot]=em_alg_function(x,m,s,Pa,e_min)
% EM algorithm for estimating the parameters of a mixture of normal
% distributions, with diagonal covariance matrices.
% WARNING: IT ONLY RUNS FOR THE CASE WHERE THE COVARIANCE MATRICES
% ARE OF THE FORM sigma^2*I. IN ADDITION, IF sigma_i^2=0 FOR SOME
% DISTRIBUTION AT AN ITERATION, IT IS ARBITRARILY SET EQUAL TO 0.001.
%
% INPUT ARGUMENTS:
%   x:      lxN matrix, each column of which is a feature vector.
%   m:      lxJ matrix, whos j-th column is the initial
%           estimate for the mean of the j-th distribution.
%   s:      1xJ vector, whose j-th element is the variance
%           for the j-th distribution.
%   Pa:     J-dimensional vector, whose j-th element is the initial
%           estimate of the a priori probability of the j-th distribution.
%   e_min:  threshold used in the termination condition of the EM
%           algorithm.
%
% OUTPUT ARGUMENTS:
%   m:      it has the same structure with input argument m and contains
%           the final estimates of the means of the normal distributions.
%   s:      it has the same structure with input argument s and contains
%           the final estimates of the variances of the normal
%           distributions.
%   Pa:     J-dimensional vector, whose j-th element is the final estimate
%           of the a priori probability of the j-th distribution.
%   iter:   the number of iterations required for the convergence of the
%           EM algorithm.
%   Q_tot:  vector containing the likelihood value at each iteration.
%   e_tot:  vector containing the error value at each itertion.
%
% (c) 2010 S. Theodoridis, A. Pikrakis, K. Koutroumbas, D. Cavouras
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 8
    maxiter = Inf;
end

p = length(x);
J = length(mu);

e = e_min + 1;

Q_tot = [];
e_tot = [];

global CONDPROB
global THETA

iter = 0;
while e > e_min && iter < maxiter
    iter = iter + 1;
    e
    
    P_old = Pa;
    mu_old = mu;
    kappa_old = kappa;
    
    % Determine P(j|x_k; theta(t))
    for k = 1:p
        for j = 1:J
            temp(j) = vmpdf(x(k),mu(j),kappa(j));
        end
        P_tot = temp * Pa';
        for j = 1:J
            P(j,k) = temp(j) * Pa(j) / P_tot;
        end
    end
    
    % Determine the log-likelihood
%     Q=0;
%     for k = 1:p
%         L = 0;
%         for j = 1:J
%             L0 = Pa(j) / (2 * pi * besseli(0,kappa(j))) * exp(kappa(j)*cos(x(k)-mu(j)));
%             L = L + L0;
%         end
%         Q = Q + P(j,k) * log(L);
%     end
%     Q_tot = [Q_tot Q];
    CONDPROB = P;
    THETA = x;
    
    [prms,f] = fminsearchbnd(@loglh,[mu kappa Pa(1:end-1)],bl,bu);
    Q_tot = [Q_tot -f];
    
    % Estimated parameters
    mu = prms(1:J);
    kappa = prms(J+1:2*J);
    Pa = prms(2*J+1:end);
    Pa(end+1) = 1 - sum(Pa);
    
%     % Determine mu and kappa
%     for j = 1:J
%         a = 0;
%         for k = 1:p
%             a = a + P(j,k) / sum(P(j,:)) * exp(i*x(k));
%         end
%         mu(j) = angle(a);
%         mu(mu<0) = mu(mu<0) + 2 * pi
%         kappa(j) = A1inv(abs(a))
%     end
%     
%     % Determine the a priori probabilities
%     for j = 1:J
%         a = 0;
%         for k = 1:p
%             a = a + P(j,k);
%         end
%         Pa(j) = a / p;
%     end
    
    e = sum(abs(Pa-P_old)) + sum(abs(mu-mu_old)) + sum(abs(kappa-kappa_old));    
    e_tot = [e_tot e];
end

% -------------------------------------------------------------------------
function r = loglh(xc)
%Log-likelihood

% Define parameters
global CONDPROB
P = CONDPROB;
global THETA
x = THETA;
k = (length(xc) + 1) / 3;
mu = xc(1:k);
kappa = xc(k+1:2*k);
Pa = xc(2*k+1:end);
Pa(end+1) = 1 - sum(Pa);
p = length(x);
J = length(mu);

Q=0;
% for k = 1:p
%     L = 0;
%     for j = 1:J
%         L0 = Pa(j) / (2 * pi * besseli(0,kappa(j))) * exp(kappa(j)*cos(x(k)-mu(j)));
%         L = L + L0;
%     end
%     Q = Q + P(j,k) * log(L);
% end
% Q2 = Q;
% Q = 0;
% for k = 1:p
%     L0 = Pa(1:J) ./ (2 * pi * besseli(0,kappa(1:J))) .* exp(kappa(1:J).*cos(x(k)-mu(1:J)));
%     L = sum(L0);
%     Q = Q + P(j,k) * log(L);
% end
Q=0;
for k=1:p
    for j=1:J
        Q=Q+P(j,k)*(-log(2*pi*besseli(0,kappa(j))) + cos(x(k)-mu(j))*kappa(j) + log(Pa(j)));
    end
end
% Q_tot = [Q_tot Q];
r = -Q;