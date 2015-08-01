function [pvalue AIC BIC] = vmcomponents_EM(theta,degrad)
%VMCOMPONENTS_EM   Estimates the number of von Mises components.
%   [P AIC BIC] = VMCOMPONENTS_EM(THETA,DEGRAD) performs statistical tests
%   addressing the question whether the number of modes in the circular
%   distribution of the sample THETA is (i) k or (ii) more than k, from k =
%   1 to 5. The second input argument (DEGRAD) should contain the
%   information whether the sample is measured in degress ('deg') or
%   radians ('rad'). P values are returned along with AIC and BIC
%   information criteria for model selection.
%
%   VMCOMPONENTS_EM applies an Expectation Maximization algorithm (see
%   WATSONKFIT_EM and EM_ALG_FUNCTION_VM).
%
%   Reference: Fisher NI (1993) Statistical analysis of circular data,
%   Cambridge University Press, Cambridge. pp. 100-102.
%
%   See also WATSONKFIT and WATSONKFIT_EM.

% Input argument check
error(nargchk(2,2,nargin))
if isequal(degrad,'deg')    % convert to radian
    theta = theta / 180 * pi;
end

% Directories
global DATAPATH
inpdir = [DATAPATH 'Czurko2\VMfit\'];
resdir = [DATAPATH 'Czurko2\VMfit\bootstrap2\'];

% Estimate the number of von Mises components
pvalue = zeros(1,5);
AIC = zeros(1,5);
BIC = zeros(1,5);
for k = 1:5
    [pvalue(k) AIC(k) BIC(k)] = main(theta,k,inpdir,resdir);
    save([resdir 'pvalue.mat'],'pvalue','AIC','BIC')
end



% -------------------------------------------------------------------------
function [pvalue AIC BIC] = main(theta,k,inpdir,resdir)

% Load parameter estimations (Step 1.)
n = length(theta);
switch k
    case 1
        [ftm, mu, mvl] = mvlmn(theta,'rad');
        kappa = A1inv(mvl);
        p = [];
    case 2
        ff = [inpdir 'estEM2c'];
        load(ff)
    case 3
        ff = [inpdir 'estEM3a'];
        load(ff)
    case 4
        ff = [inpdir 'estEM4b'];
        load(ff)
    case 5
        ff = [inpdir 'estEM5a'];
        load(ff)
end
p(end+1) = 1 - sum(p);
F = cell(1,k);    % mixed von Mises cdf
cF = zeros(1,n);
for t = 1:k
    F{t} = vmcdf(theta',mu(t),kappa(t),'rad');
    cF = cF + p(t) * F{t};
end
L = 0;  % log-likelihhod function
for j = 1:k
    L = L + p(j)/(2*pi*besseli(0,kappa(j)))*exp(kappa(j)*cos(theta-mu(j)));
end
Q = sum(log(L));
AIC = 2 * (3 * k - 1) - 2 * Q;
BIC = - 2 * Q + (3 * k - 1) * log(n);

% Goodness-of-fit (Step 2.)
z = sort(cF,'ascend');
T0 = (z - (2 * (1:n) - 1) / (2 * n)) .^ 2;
T = sum(T0) - n * (mean(z) - 0.5) ^ 2 + 1 / (12 * n);   % U^2 (4.35)

% Parametric resample (Step. 3.)
B = 100;
T_star = zeros(1,B);
for bts = 1:B
    disp(['bts = ' num2str(bts)])
    Phat0 = 0;
    Phat = cumsum(p);   % (3.1)
    u = rand(1,n);      % U(0,1) (3.2)
    theta_star = zeros(1,n);
    for t = 1:n         % (3.3)
        j = find([Phat0 Phat]<=u(t),1,'first');
        theta_star(t) = b_rvm(1,mu(j),kappa(j));
    end

% Parameter estimation for the bootstrap sample (Step 4.)
    if k == 1
        [ftm_star, mu_star, mvl_star] = mvlmn(theta_star,'rad');
        kappa_star = A1inv(mvl_star);
        p_star = [];
    else
        [param,err] = watsonkfit_EM(theta_star,k);     % parameter estimation
        mu_star = param{1};
        kappa_star = param{2};
        p_star = param{3};
    end
    p_star(end+1) = 1 - sum(p_star);

    F_star = cell(1,k);    % mixed von Mises cdf
    cF_star = zeros(1,n);
    for t = 1:k
        F_star{t} = vmcdf(theta_star,mu_star(t),kappa_star(t),'rad');
        cF_star = cF_star + p_star(t) * F_star{t};
    end

    z = sort(cF_star,'ascend');  % goodness-of-fit
    T0 = (z - (2 * (1:n) - 1) / (2 * n)) .^ 2;
    T_star(bts) = sum(T0) - n * (mean(z) - 0.5) ^ 2 + 1 / (12 * n);
    save([resdir 'Tstar_' num2str(k) '.mat'],'T_star')
end

% Test k modes against more than k modes
pvalue = (1 / B) * sum(T<=T_star);