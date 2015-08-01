function [param,err] = b_watsontwofit2(theta,pl)
%WATSONTWOFIT2   Fits mixture of two von Mises distributions.
%   [PARAM, ERR] = WATSONTWOFIT2(THETA) fits a mixture of two von Mises
%   distributions on the circular data set THETA. Note that THETA should
%   contain radian values! PARAM is a 1-by-5 array containing mu1, mu2,
%   kappa1, kappa2 - mean and concentration parameter of the two von Mises
%   distributions - and p mixing ratio (0 < p < 1). ERR stands for error
%   function value at the optimum of minimization.
%
%   [PARAM, ERR] = WATSONTWOFIT2(THETA,PL) plots result, if PL = 1.
%
%   As compared to WATSONTWOFIT, it uses a slightly different minimization
%   algorithm: it minimizes parameters first at fixed p = 0.5, then it
%   minimizes p while fixing parameters to the found optimum.
%
%   See also WATSONTWOFIT, WATSON, WATSONTWO, RVM and RMIXEDVM.

% Input argument check
error(nargchk(1,2,nargin))
if nargin == 1
    pl = 0;
end
theta = mod(theta,2*pi);

% Declare 'THETA' as global variable in order to pass to 'funct3'
global THETA
THETA = theta;

% Estimate least square minimum
bl = [0 0 0.3 0.3];       % lower bound
bu = [2*pi 2*pi 10 10]; % upper bound
[x,f] = fminsearchbnd('funct3',[1,2,0.5,0.5],bl,bu);

% Estimated parameters
mu1 = x(1);
mu2 = x(2);
kappa1 = x(3);
kappa2 = x(4);

% Pass parameters
global PARAM
PARAM(1) = mu1;
PARAM(2) = mu2;
PARAM(3) = kappa1;
PARAM(4) = kappa2;

% Estimate 'p'
[p,ff] = fminsearchbnd('funct4',[0.5],0,1);

% Clear global variables
clear global THETA
clear global PARAM

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