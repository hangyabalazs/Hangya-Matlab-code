function [v ft t ITIs] = bimodal_subjective_hazard(Phi,dt,NumTrials)
%BIMODAL_SUBJECTIVE_HAZARD   Bimodal subjective hazard rate model.
%   [HWT FT T] = BIMODAL_HAZARD_MODEL(PHI,DT) calculates the subjective
%   hazard rate (HWT) and the failure density function (FT) based on the
%   bimodal anticipation function (T, corresponding time vector). Time step
%   is determined by DT. The bimodal anticipation function is a mixture of
%   two Gaussians (means, 0.3 and 2s; SD, 0.15s) and a uniform random
%   variable constrained between 0.1 and 3 seconds (mixing proportions:
%   0.35, 0.35, 0.3). Subjective hazard rate is calculated from the hazard
%   rate by convolving it with a Gaussian pdf with variability linearly
%   increasing in time (Phi, slope parameter of this increase).
%
%   [HWT FT T ITIS] = BIMODAL_SUBJECTIVE_HAZARD(PHI,DT,NUMTRIALS) returns a
%   random sample of failure times (ITIS; NUMTRIALS, sample size).
%
%   Reference:
%   Janssen P, Shadlen MN (2005) A representation of the hazard rate of
%   elapsed time in macaque area LIP. Nat Neurosci 8(2):234-41.
%
%   See alse BIMODALITI_HAZARD_MODEL.

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
if nargin < 2
    dt = 0.01;     % resolution
end
t = ITIMin:dt:ITIMax;   % time

% Simulate
if nargout == 4
    if nargin < 3
        NumTrials = 1000;
    end
    
    % Sample the failure distribution
    ITIs1 = random('Normal',mng1,sdg,1,NumTrials);
    while any(ITIs1>ITIMax) || any(ITIs1<ITIMin)
        inx = ITIs1 > ITIMax  | ITIs1 < ITIMin;
        ITIs1(inx) = random('Normal',mng1,sdg,1,sum(inx));
    end
    ITIs2 = random('Normal',mng2,sdg,1,NumTrials);
    while any(ITIs2>ITIMax) || any(ITIs2<ITIMin)
        inx = ITIs2 > ITIMax  | ITIs2 < ITIMin;
        ITIs2(inx) = random('Normal',mng2,sdg,1,sum(inx));
    end
    ITIs3 = random('Uniform',ITIMin,ITIMax,1,NumTrials);
    prr = rand(1,NumTrials);
    rr = zeros(3,NumTrials);
    rr(1,prr<pmx1) = 1;
    rr(2,prr>=pmx1&prr<(pmx1+pmx2)) = 1;
    rr(3,prr>=(pmx1+pmx2)) = 1;
    ITIs = rr(1,:) .* ITIs1 + rr(2,:) .* ITIs2 + rr(3,:) .* ITIs3;
end

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
pret = fliplr(t(1)-dt:-dt:dt);
postt = t(end)+dt:dt:10*t(end);
t2 = [pret t postt];   % we need a longer time range, to get a better estimate of Fwt
lenpt = length(pret);
tinx = lenpt+1:lenpt+length(t);   % indices for the original times
for tm = t2
    fg = ft .* exp(-(tau-tm).^2/(2*Phi^2*tm^2));
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

% Values corresponding to the original times
v = hwt(tinx);