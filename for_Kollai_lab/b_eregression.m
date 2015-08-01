function b_eregression(data_control,data_patient)
%EREGRESSION   Forward Stepwise Regression for 'Ellasticity' data.
%   EREGRESSION(D1,D2) calculates Forward Stepwise Regression for control (D1) and 
%   patient (D2) 'Ellasticity' data file. The regressors are systolic blood pressure,
%   diastolic blood pressure, age, heart rate and weight. The dependent variable is
%   complience, distensibility coefficient, stiffness index or incremental ellastic
%   modulus. The function saves the regression results in excel files and the residuals
%   in mat files. You are able to modify the destination directory through editing 
%   the program code.
%
%   See also EIMPORT, EANOCOVA, ECOMPARE, ECOMPARE2, ECOMPARE3 and EREGRESSION2.

% Input argument check
error(nargchk(2,2,nargin))

% Define output directory and filenames
outdir = 'F:\balazs\MersichBea2\';
xlsfile_control = [outdir 'regression_control.xls'];
xlsfile_patient = [outdir 'regression_patient.xls'];
matfile_control = [outdir 'residuals_control.mat'];
matfile_patient = [outdir 'residuals_patient.mat'];

% Getting the regressors
sbp_control = data_control(:,9);    % systolic blood pressure
sbp_patient = data_patient(:,9);
dbp_control = data_control(:,10);   % diastolic blood pressure
dbp_patient = data_patient(:,10);
age_control = data_control(:,4);    % age
age_patient = data_patient(:,4);
hr_control = data_control(:,16);    % heart rate
hr_patient = data_patient(:,16);
weight_control = data_control(:,7); % weight
weight_patient = data_patient(:,7);
comp_control = data_control(:,22);  % complience
comp_patient = data_patient(:,22);
dc_control = data_control(:,23);    % distensibility coefficient
dc_patient = data_patient(:,23);
stif_control = data_control(:,24);  % stiffness index
stif_patient = data_patient(:,24);
einc_control = data_control(:,27);  % incremental ellastic modulus
einc_patient = data_patient(:,27);

X_control = [sbp_control dbp_control age_control hr_control weight_control];
X_patient = [sbp_patient dbp_patient age_patient hr_patient weight_patient];

% Forward Stepwise Regression
res_control_comp = fsr(X_control,comp_control,'Control',xlsfile_control,'comp');
res_control_dc = fsr(X_control,dc_control,'Control',xlsfile_control,'dc');
res_control_stif = fsr(X_control,stif_control,'Control',xlsfile_control,'stif');
res_control_einc = fsr(X_control,einc_control,'Control',xlsfile_control,'einc');
res_patient_comp = fsr(X_patient,comp_patient,'Patient',xlsfile_patient,'comp');
res_patient_dc = fsr(X_patient,dc_patient,'Patient',xlsfile_patient,'dc');
res_patient_stif = fsr(X_patient,stif_patient,'Patient',xlsfile_patient,'stif');
res_patient_einc = fsr(X_patient,einc_patient,'Patient',xlsfile_patient,'einc');

% Saving the residuals
save(matfile_control,'res_control_comp','res_control_dc','res_control_stif','res_control_einc')
save(matfile_patient,'res_patient_comp','res_patient_dc','res_patient_stif','res_patient_einc')

% -----------------------------------------------------------------------------------------
function r = fsr(X,y,corp,filename,depvarname)
%FSR   Forward Stepwise Regression.
%   Input arguments:
%       X: regressors
%       y: data
%       corp: 'Control' or 'Patient'
%       filename: name of output file
%       depvarname: name of dependent variable
%   Output argument:
%       res: array of residuals

% Regression
stepwise(X,y)
[bcoeff,se,pval,inmodel,stats,nextstep,history] = stepwisefit(X,y,'display','off');
hin = history.in;
hin(end+1,:) = ones(1,size(hin,2));
ord = zeros(5,1);
for i = 1:size(history.in,1)
    fnd = find(hin(i,:));
    ord(fnd) = i;
    hin(i+1,:) = hin(i+1,:) - history.in(i,:);
end
tstat = stats.TSTAT;
titlestr = {corp depvarname};
headers = {'bcoeff' 'stand. err.' 't-stat' 'p' 'order'};
regressors = {'sbp' 'dbp' 'age' 'hr' 'weight'}';
warning off MATLAB:xlswrite:AddSheet
xlswrite(filename,titlestr,depvarname);     % generating the excel file
xlswrite(filename,headers,depvarname,'B2');
xlswrite(filename,regressors,depvarname,'A3');
xlswrite(filename,[bcoeff se tstat pval ord],depvarname,'B3');

% Calculating residuals
X1 = [ones(size(X,1),1) X];
[b,bint,r,rint,regstats] = regress(y,X1);