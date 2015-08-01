function [param err] = fitHazardRate(initparam,iti,rt,hazardmodelfun)
%FITHAZARDRATE   Fit hazard rate model on behavioral data.
%   [PARAM,ERR] = FITHAZARDRATE(INITPARAM,ITI,RT,HAZARDMODELFUN) fits
%   subjective hazard rate on sample foreperiod (ITI) and reaction time
%   (RT) using the model HAZARDMODELFUN (function handle). Initial model
%   parameters shoud be passed in INITPARAM and optimal parameters are
%   returned in PARAM. The fitting is perfomred by the log-likelihood
%   method and the error of the fit (negative log-likelihood) is returned
%   in ERR. The optimization is implemented using FMINSEARCH (Nelder-Mead
%   simplex algorithm).
%
%   See alse ERRORHAZARDRATE, BIMODALITI_FIT and BIMODAL_HAZARD_MODEL.

% Optimize using 'fminsearch'
[param err] = fminsearchcall(@(p)errorHazardRate(p,iti,rt,hazardmodelfun),...
    initparam,optimset('MaxFunEvals',150));

% -------------------------------------------------------------------------
function [outparam err] = fminsearchcall(errorfun,initparam,opt)

% Transform parameter input to a vector
newinitparam = [];
fld = fieldnames(initparam);
lenfld = length(fld);
for k = 1:lenfld
    newinitparam = [newinitparam initparam.(fld{k})]; %#ok<AGROW>
end

% Call 'fminsearch'
[param err] = fminsearch(errorfun,newinitparam,opt);

% Transform output parameter back to structure
for k = 1:lenfld
    outparam.(fld{k}) = param(k);
end