function b_eregression2(data_control,data_patient)
%EREGRESSION2   Linear Regression for 'Ellasticity' data.
%   EREGRESSION2(D1,D2) calculates Linear Regression for control (D1) and patient 
%   (D2) 'Ellasticity' data file. The only regressor is systolic blood pressure.
%   The dependent variable is complience, distensibility coefficient or incremental 
%   ellastic modulus. The function saves the residuals in mat files. You are able
%   to modify the destination directory through editing the program code.
%
%   See also EIMPORT, EANOCOVA, ECOMPARE, ECOMPARE2 and ECOMPARE3.

% Input argument check
error(nargchk(2,2,nargin))

% Define output directory and filenames
outdir = 'F:\balazs\MersichBea3\';
matfile_control = [outdir 'residuals_control.mat'];
matfile_patient = [outdir 'residuals_patient.mat'];

% Getting the regressors
sbp_control = data_control(:,9);    % systolic blood pressure
sbp_patient = data_patient(:,9);

comp_control = data_control(:,22);  % complience
comp_patient = data_patient(:,22);
dc_control = data_control(:,23);    % distensibility coefficient
dc_patient = data_patient(:,23);
einc_control = data_control(:,27);  % incremental ellastic modulus
einc_patient = data_patient(:,27);

% Linear Regression
res_control_comp = lr(sbp_control,comp_control,'Control','comp');
res_control_dc = lr(sbp_control,dc_control,'Control','dc');
res_control_einc = lr(sbp_control,einc_control,'Control','einc');
res_patient_comp = lr(sbp_patient,comp_patient,'Patient','comp');
res_patient_dc = lr(sbp_patient,dc_patient,'Patient','dc');
res_patient_einc = lr(sbp_patient,einc_patient,'Patient','einc');

% Saving the residuals
save(matfile_control,'res_control_comp','res_control_dc','res_control_einc')
save(matfile_patient,'res_patient_comp','res_patient_dc','res_patient_einc')

% -----------------------------------------------------------------------------------------
function r = lr(X,y,corp,depvarname)
%FSR   Linear Regression.
%   Input arguments:
%       X: regressors
%       y: data
%       corp: 'Control' or 'Patient'
%       depvarname: name of dependent variable
%   Output argument:
%       res: array of residuals

% Regression and residuals
X1 = [ones(size(X,1),1) X];
[b,bint,r,rint,regstats] = regress(y,X1);