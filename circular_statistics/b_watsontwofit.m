function [param,err] = b_watsontwofit(theta,pl)
%WATSONTWOFIT   Fits mixture of two von Mises distributions.
%   [PARAM, ERR] = WATSONTWOFIT(THETA) fits a mixture of two von Mises
%   distributions on the circular data set THETA. Note that THETA should
%   contain radian values! PARAM is a 1-by-5 array containing mu1, mu2,
%   kappa1, kappa2 - mean and concentration parameter of the two von Mises
%   distributions - and p mixing ratio (0 < p < 1). ERR stands for error
%   function value at the optimum of minimization.
%
%   [PARAM, ERR] = WATSONTWOFIT(THETA,PL) plots result, if PL = 1 or if PL
%   is undefined.
%
%   See also WATSON, WATSONTWO, RVM and RMIXEDVM.

% Input argument check
error(nargchk(1,2,nargin))
if nargin == 1
    pl = 1;
end
theta = mod(theta,2*pi);

% Declare 'THETA' as global variable in order to pass to 'funct2'
global THETA
THETA = theta;

% Estimate least square minimum
bl = [0 0 0.3 0.3 0];       % lower bound
bu = [2*pi 2*pi 10 10 1];   % upper bound

kappa1 = 5;     % initial estimation
kappa2 = 5;
p = 0.5;
next = 1;
for mu1 = 0:0.1:2*pi
    for mu2 = mu1:0.1:2*pi
        err.f(next) = funct2([mu1 mu2 kappa1 kappa2 p]);
        err.mu1(next) = mu1;
        err.mu2(next) = mu2;        
        next = next + 1;
    end
end
inx = find(err.f==min(err.f));
mu1 = err.mu1(inx);
mu2 = err.mu2(inx);

[x,f] = fminsearchbnd('funct2',[mu1,mu2,kappa1,kappa2,p],bl,bu);
if f > 0.02
    disp(['Estimation error: ' f]);
end

% Estimated parameters
mu1 = x(1);
mu2 = x(2);
kappa1 = x(3);
kappa2 = x(4);
p = x(5);

% Clear global 'THETA'
clear global THETA

% Plot von Mises
if pl
%     figure
    x = linspace(0,2*pi,length(theta));
    dvm = p / (2 * pi * besseli(0,kappa1,1)) * (exp(cos(x-mu1)-1)).^kappa1 +...
        (1 - p) / (2 * pi * besseli(0,kappa2,1)) * (exp(cos(x-mu2)-1)).^kappa2;
    hold on
    bno = fix(exp(0.626+0.4*log(length(theta)-1)));   % number of bins
    [thb,nb] = hist(theta,bno);
    thb_norm = thb / mean(thb) * mean(dvm);
    bar([nb nb+2*pi],[thb_norm thb_norm]);
    plot([x x+2*pi],[dvm dvm],'r','LineWidth',2)
end

% Output arguments
param(1) = mu1;
param(2) = mu2;
param(3) = kappa1;
param(4) = kappa2;
param(5) = p;
err = f;