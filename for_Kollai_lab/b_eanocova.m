function b_eanocova(data_control,data_patient)
%EANOCOVA   Analysis of covariance for 'Ellasticity' data.
%   EANOCOVA(D1,D2) runs AOCTOOL for D1 and D2 choosing systolic blood
%   pressure as independent and complience, distensibility coefficient,
%   stiffness index or incremental ellasticity modulus as dependent
%   variable.
%
%   See also EIMPORT, EREGRESSION, EREGRESSION2, ECOMPARE, ECOMPARE2,
%   ECOMPARE3, and AOCTOOL.

% Input argument check
error(nargchk(2,2,nargin))

% Define output directory and filenames
outdir = 'F:\balazs\MersichBea2\';

% Getting the regressors
sbp_control = data_control(:,9);    % systolic blood pressure
sbp_patient = data_patient(:,9);
sbp = [sbp_control; sbp_patient];
dbp_control = data_control(:,10);   % diastolic blood pressure
dbp_patient = data_patient(:,10);
dbp = [dbp_control; dbp_patient];
age_control = data_control(:,4);    % age
age_patient = data_patient(:,4);
age = [age_control; age_patient];
hr_control = data_control(:,16);    % heart rate
hr_patient = data_patient(:,16);
hr = [hr_control; hr_patient];
weight_control = data_control(:,7); % weight
weight_patient = data_patient(:,7);
weight = [weight_control; weight_patient];
comp_control = data_control(:,22);  % complience
comp_patient = data_patient(:,22);
comp = [comp_control; comp_patient];
dc_control = data_control(:,23);    % distensibility coefficient
dc_patient = data_patient(:,23);
dc = [dc_control; dc_patient];
stif_control = data_control(:,24);  % stiffness index
stif_patient = data_patient(:,24);
stif = [stif_control; stif_patient];
einc_control = data_control(:,27);  % incremental ellastic modulus
einc_patient = data_patient(:,27);
einc = [einc_control; einc_patient];

group = [zeros(size(sbp_control)); ones(size(sbp_patient))];        % 0: control    1: patient

[h_comp,atab_comp,ctab_comp,stats_comp] = aoctool(sbp,comp,group);
[h_dc,atab_dc,ctab_dc,stats_dc] = aoctool(sbp,dc,group);
[h_einc,atab_einc,ctab_einc,stats_einc] = aoctool(sbp,einc,group);