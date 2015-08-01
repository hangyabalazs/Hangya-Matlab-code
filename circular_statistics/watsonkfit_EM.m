function [param,err,Q_tot,e_tot] = watsonkfit_EM(theta,k,pl)
%WATSONKFIT_EM   Fits mixture of k von Mises distributions with Expectation Maximization.
%   [PARAM, ERR] = WATSONKFIT_EM(THETA,K) fits a mixture of K von Mises
%   distributions on the circular data set THETA. Note that THETA should
%   contain radian values! PARAM returns mean and concentration parameter
%   values of the K von Mises distributions as well as the  mixing ratios.
%   ERR stands for error function value in the optimum of minimization.
%
%   WATSONKFIT_EM applies an Expectation Maximization aldorithm implemented
%   for circular von Mises distributions.
%
%   [PARAM, ERR, Q, E] = WATSONKFIT_EM(THETA,K) returns log-likelihood (Q)
%   and error (E) values for all iterations.
%
%   [PARAM, ERR] = WATSONKFIT_EM(THETA,PL) plots result, if PL = 1 or if PL
%   is undefined.
%
%   See also WATSONTWOFIT, WATSONKFIT and EM_ALG_FUNCTION_VM.

% Input argument check
error(nargchk(2,3,nargin))
if nargin == 2
    pl = 0;
end
theta = mod(theta,2*pi);

% Declare 'THETA' as global variable in order to pass to 'funct2'
global THETA
THETA = theta;

% Estimate least square minimum
bl = [zeros(1,k) ones(1,k)*0.3 zeros(1,k-1)];       % lower bound
bu = [ones(1,k)*2*pi ones(1,k)*25 ones(1,k-1)];   % upper bound

kappa = ones(1,k) * 10;     % initial estimation
% p = ones(1,k-1) / k;

next = 1;
cmb_mu = combnk(0:0.2:2*pi,k);
% cmb_kappa = combnk(8:2:16,k);
% cmb_kappa = [8 8 8 8; 8 8 8 16; 8 8 16 8; 8 16 8 8; 16 8 8 8; 8 8 16 16; 8 16 8 16; 16 8 8 16;...
%     8 16 16 8; 16 8 16 8; 16 16 8 8];
% cmb_kappa = [10 10 10 10];
cmb_kappa = kappa;
lc_mu = size(cmb_mu,1);
lc_kappa = size(cmb_kappa,1);
err = struct();
err.f = zeros(1,lc_mu*lc_kappa*10);
err.mu = cell(1,lc_mu*lc_kappa*10);
err.kappa = cell(1,lc_mu*lc_kappa*10);
for ps = 1:10
    disp(num2str(ps))
    p = rand(1,k-1) / (k -1);
    for pkp = 1:lc_kappa
        kappa = cmb_kappa(pkp,:);
        for t = 1:lc_mu
            mu = cmb_mu(t,:);
            err.f((ps-1)*lc_mu*lc_kappa+(pkp-1)*lc_mu+t) = kfunct2([{mu} {kappa} {p}]);
            err.mu{(ps-1)*lc_mu*lc_kappa+(pkp-1)*lc_mu+t} = mu;
            err.kappa{(ps-1)*lc_mu*lc_kappa+(pkp-1)*lc_mu+t} = kappa;
            next = next + 1;
        end
    end
end
inx = find(err.f==min(err.f));
mu = err.mu{inx};
kappa = err.kappa{inx};

[mu,kappa,p,iter,Q_tot,e_tot,P] = em_alg_function_VM(theta,mu,kappa,[p 1-sum(p)],bl,bu,0.02); % EM
f = e_tot(end);
if f > 0.02
    disp(['Estimation error: ' f]);
end

% Estimated parameters
p = p(1:end-1);

% Clear global 'THETA'
clear global THETA

% Plot von Mises
if pl
    figure
    x = linspace(0,2*pi,length(theta));
    q = p;
    q(end+1) = 1 - sum(q);
    pdvm = zeros(k,length(x));
    for t = 1:k
        pdvm(t,:) = q(t) .* (1 ./ (2 * pi * besseli(0,kappa(t),1)) .* (exp(cos(x-mu(t))-1)).^kappa(t));
    end
    dvm = sum(pdvm);
    hold on
    bno = fix(exp(0.626+0.4*log(length(theta)-1)));   % number of bins
    [thb,nb] = hist(theta,bno);
    thb_norm = thb / mean(thb) * mean(dvm);
    bar([nb nb+2*pi],[thb_norm thb_norm]);
    plot([x x+2*pi],[dvm dvm],'r','LineWidth',2)
    
    figure
    hold on
    edges=[(0:20:360)/180*pi];
    [thb,nb] = histc(theta,edges);
    thb_norm = thb / mean(thb) * mean(dvm);
    bar([edges(1:end-1)+10/180*pi edges(1:end-1)+10/180*pi+2*pi],[thb_norm(1:end-1)' thb_norm(1:end-1)']);
    plot([x x+2*pi],[dvm dvm],'r','LineWidth',2)
end

% Output arguments
param{1} = mu;
param{2} = kappa;
param{3} = p;
err = f;