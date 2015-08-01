function y = bimodal_hazard_model2(param,x)
%BIMODAL_HAZARD_MODEL2   Bimodal subjective hazard rate model.
%   RT = BIMODAL_HAZARD_MODEL2(PARAM,ITI) calculates model results (RT,
%   reation time) at given points (ITI, foreperiod) based on the following
%   subjective hazard rate model:
%       RT(t) = we + wr * exp(x) + wb * Ab(t)
%   where 'we' is a constant parameter, 'wr' and 'wb' are the weights of
%   the exponential function and the bimodal anticipation function denoted
%   by 'Ab', respectively. Model parameters should be passed to the
%   function in the PARAM structure. The bimodal anticipation function (Ab)
%   is a mixture of two Gaussians (means, 0.3 and 2s; SD, 0.15s) and a
%   uniform random variable constrained between 0.1 and 3 seconds (mixing
%   proportions: 0.35, 0.35, 0.3). Subjective hazard rate is calculated
%   from the hazard rate by convolving it with a Gaussian pdf with
%   variability linearly increasing in time (Phi, slope parameter of this
%   increase).
%
%   Reference:
%   Janssen P, Shadlen MN (2005) A representation of the hazard rate of
%   elapsed time in macaque area LIP. Nat Neurosci 8(2):234-41.
%
%   See alse FITHAZARDRATE, ERRORHAZARDRATE and BIMODALITI_FIT.

% Handle both vector and structure input
if ~isstruct(param)
    st.we = param(1);
    st.wr = param(2);
    st.wb = param(3);
    st.phi = param(4);
    st.tau = param(5);
    param = st;
end

% Call the model
y = param.we + param.wr * exp(x-param.tau) +...
    param.wb * subjective_hazard_rate(param.phi,x-param.tau);

% -------------------------------------------------------------------------
function v = subjective_hazard_rate(Phi,iti)

% Parameters of the foreperiod distribution
ITIMin = 0.1;   % parameters of the uniform distribution
ITIMax = 3;
mng1 = 0.3;   % parameters for the Gaussians
mng2 = 2;
sdg = 0.15;
pmx1 = 0.35;   % mixing probabilities
pmx2 = 0.35;
pmx3 = 1 - pmx1 - pmx2;

% Time
dt = 0.01;     % resolution
t = ITIMin:dt:ITIMax;   % time

% Failure density function
ft = pmx3 / (ITIMax - ITIMin) + ...
    pmx1 * normpdf(t,mng1,sdg) + ...
    pmx2 * normpdf(t,mng2,sdg);
% figure
% plot(t,ft)
ft = ft / sum(ft);  % for estimating cdf, we need a histogram estimation

% Hazard rate
Ft = cumsum(ft);
Rt = 1 - Ft;
ht = (Rt(1:end-1) - Rt(2:end)) ./ (dt * Rt(1:end-1));
% figure
% plot(t(1:end-1),ht)

% Subjective hazard rate:
% Convolution of the failure density function with a time-dependent Gaussian
intfg = zeros(size(t));
tau = t;
cntr = 1;   % counter
t2 = [fliplr(t(1)-dt:-dt:dt) t t(end)+dt:dt:3*t(end)];   % we need a longer time range, to get a better estimate of Fwt
for tm = t2
    fg = ft .* exp(-(tau-tm).^2/(2*(Phi*tm).^2));
    intfg(cntr) = sum(fg*dt);
    cntr = cntr + 1;
end
fwt = 1 ./ (Phi * t2 * sqrt(2*pi)) .* intfg;
fwt = fwt / sum(fwt);

% Subjective hazard rate
Fwt = cumsum(fwt);
Rwt = 1 - Fwt;
hwt = (Rwt(1:end-1) - Rwt(2:end)) ./ (dt * Rwt(1:end-1));
% figure
% plot(t2(1:end-1),hwt)

% Values at requested times
v = linterp(t2(1:end-1),hwt,iti);