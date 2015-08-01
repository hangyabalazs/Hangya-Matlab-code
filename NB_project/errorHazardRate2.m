function R = errorHazardRate2(param,iti,rt,hazardmodelfun)
%ERRORHAZARDRATE2   Least square error for subjective hazard rate model.
%   R = ERRORHAZARDRATE2(PARAM,ITI,RT,HAZARDMODELFUN) calculates least
%   square error (R) for the sample reaction time (RT) and corresponding
%   foreperiod (ITI) data based on hazard rate model HAZARDMODELFUN
%   (function handle).Parameters should be passed to the function in the
%   PARAM structure.
%
%   See alse FITHAZARDRATE2, BIMODAL_HAZARD_MODEL and BIMODALITI_FIT.

% Model mean values
y_hat = feval(hazardmodelfun,param,iti);

% Least square error
R = nansum((y_hat-rt).^2);
disp(R)