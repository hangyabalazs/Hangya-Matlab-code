function b_ecompare3
%ECOMPARE3   Comparation of 'Ellasticity' data regression residuals.
%   ECOMPARE3 calculates F-test for residuals saved by EREGRESSION. The
%   function saves the results (test statistics values and p values) in a
%   text file, which you are able to modify through editing the program
%   code.
%
%   See also EIMPORT, EREGRESSION, EREGRESSION2, EANOCOVA, ECOMPARE and 
%   ECOMPARE2.

% Input argument check
error(nargchk(0,0,nargin))

% Define output filename
outdir = 'F:\balazs\MersichBea3\';
outfile = [outdir 'Fstat.txt'];

% Load residuals
inpdir = 'F:\balazs\MersichBea2\';
matfile_control = [outdir 'residuals_control.mat'];
matfile_patient = [outdir 'residuals_patient.mat'];
load(matfile_control)
load(matfile_patient)

% F-test
[Fstat_comp,Fp_comp] = Ftest(res_control_comp,res_patient_comp);
[Fstat_dc,Fp_dc] = Ftest(res_control_dc,res_patient_dc);
[Fstat_einc,Fp_einc] = Ftest(res_control_einc,res_patient_einc);

% Save
fid = fopen(outfile,'w');
fprintf(fid,'%s\n','comp');
fprintf(fid,'%s %f %s %f\n','F: ',Fstat_comp,'              p: ',Fp_comp);
fprintf(fid,'\n');

fprintf(fid,'%s\n','dc');
fprintf(fid,'%s %f %s %f\n','F: ',Fstat_dc,'              p: ',Fp_dc);
fprintf(fid,'\n');

fprintf(fid,'%s\n','einc');
fprintf(fid,'%s %f %s %f\n','F: ',Fstat_einc,'              p: ',Fp_einc);
fprintf(fid,'\n');
fclose(fid);

% -------------------------------------------------------------------------
function [stat,p] = Ftest(d1,d2)
%FTEST   F-test.
%   [STAT,P] = FTEST(D1,D2) returns the F-stat. (STAT) and the
%   corresponding p (P) value for D1 and D2.

s1 = max(var(d1),var(d2));
s2 = min(var(d1),var(d2));
stat = s1 / s2;
v1 = length(d1) - 1;
v2 = length(d2) - 1;
if isequal(s1,var(d1))
    p = 2 * (1 - fcdf(stat,v1,v2));
elseif isequal(s1,var(d2))
    p = 2 * (1 - fcdf(stat,v2,v1));
end