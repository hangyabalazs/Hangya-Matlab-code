function L = errorHazardRate(param,iti,rt,hazardmodelfun)
%ERRORHAZARDRATE   Negative log-likelihood for subjective hazard rate model.
%   L = ERRORHAZARDRATE(PARAM,ITI,RT,HAZARDMODELFUN) calculates negative
%   log-likelihood (L) for the sample reaction time (RT) and corresponding
%   foreperiod (ITI) data based on hazard rate model HAZARDMODELFUN
%   (function handle). The model distribution is constructed as follows:
%       RT(t) = HAZARDMODELFUN(t) + E
%   where E is a Gaussian error variable with zero mean and standard
%   deviation of the sample reation time. Parameters should be passed to
%   the function in the PARAM structure.
%
%   See alse FITHAZARDRATE, BIMODAL_HAZARD_MODEL and BIMODALITI_FIT.

% Model mean values
y_hat = feval(hazardmodelfun,param,iti);

% Probability of the sample
y_prob = normpdf(rt,y_hat,nanstd(rt));     % assuming Gaussian distribution around the mean, with variance derived from the sample

% Negative log-likelihood
L = -nansum(log(y_prob));
disp(L)