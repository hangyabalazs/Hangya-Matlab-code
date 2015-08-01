function b_ecompare
%ECOMPARE   Comparation of 'Ellasticity' data regression residuals.
%   ECOMPARE calculates F-test and continuous homogenity test
%   (Kolmogorov-Smirnov-test) for residuals saved by EREGRESSION. The
%   function saves the results (test statistics values and p values) in a
%   text file, which you are able to modify through editing the program
%   code.
%
%   See also EIMPORT, EREGRESSION, EREGRESSION2, EANOCOVA, ECOMPARE2 and 
%   ECOMPARE3.

% Input argument check
error(nargchk(0,0,nargin))

% Define output filename
outdir = 'F:\balazs\MersichBea\';
outfile = [outdir 'FandKSstat.txt'];

% Load residuals
inpdir = 'F:\balazs\MersichBea\';
matfile_control = [outdir 'residuals_control.mat'];
matfile_patient = [outdir 'residuals_patient.mat'];
load(matfile_control)
load(matfile_patient)

% F-test
[Fstat_comp,Fp_comp] = Ftest(res_control_comp,res_patient_comp);
[Fstat_dc,Fp_dc] = Ftest(res_control_dc,res_patient_dc);
[Fstat_stif,Fp_stif] = Ftest(res_control_stif,res_patient_stif);
[Fstat_einc,Fp_einc] = Ftest(res_control_einc,res_patient_einc);

% Continuous homogenity test (Kolmogorov-Smirnov-test)
[KSstat_comp,KSp_comp] = KStest(res_control_comp,res_patient_comp);
[KSstat_dc,KSp_dc] = KStest(res_control_dc,res_patient_dc);
[KSstat_stif,KSp_stif] = KStest(res_control_stif,res_patient_stif);
[KSstat_einc,KSp_einc] = KStest(res_control_einc,res_patient_einc);

% Save
fid = fopen(outfile,'w');
fprintf(fid,'%s\n','comp');
fprintf(fid,'%s %f %s %f\n','F: ',Fstat_comp,'              p: ',Fp_comp);
fprintf(fid,'%s %f %s %s\n','Dn: ',KSstat_comp,'              p: ',KSp_comp);
fprintf(fid,'\n');

fprintf(fid,'%s\n','dc');
fprintf(fid,'%s %f %s %f\n','F: ',Fstat_dc,'              p: ',Fp_dc);
fprintf(fid,'%s %f %s %s\n','Dn: ',KSstat_dc,'              p: ',KSp_dc);
fprintf(fid,'\n');

fprintf(fid,'%s\n','stif');
fprintf(fid,'%s %f %s %f\n','F: ',Fstat_stif,'              p: ',Fp_stif);
fprintf(fid,'%s %f %s %s\n','Dn: ',KSstat_stif,'              p: ',KSp_stif);
fprintf(fid,'\n');

fprintf(fid,'%s\n','einc');
fprintf(fid,'%s %f %s %f\n','F: ',Fstat_einc,'              p: ',Fp_einc);
fprintf(fid,'%s %f %s %s\n','Dn: ',KSstat_einc,'              p: ',KSp_einc);
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

% -------------------------------------------------------------------------
function [Dnn,p] = KStest(d1,d2)
%KSTEST   Kolmogorov-Smirnov-test.
%   [STAT,P] = KSTEST(D1,D2) returns the Dnn stat. (STAT) and the
%   corresponding p (P) value for D1 and D2. Note, that D1 and D2 should be
%   of equally length.

n = length(d1);
if ~isequal(n,length(d2))
    error('KStest needs input arguments of equally length.')
end
Z = sort([d1; d2]');
E = zeros(size(Z));
S = zeros(size(Z));
for i = 1:length(E)
    k1 = find(d1==Z(i));
    k2 = find(d2==Z(i));
    if ~isempty(k1) & isempty(k2)
        E(i) = 1;
    elseif isempty(k1) & ~isempty(k2)
        E(i) = -1;
    else
        error('Error in KStest.')
    end
    S(i) = S(max(1,i-1)) + E(i);
end
Dnn = max(abs(S)) / n;
nDnn = Dnn * n;

if n < 4
    error('Not enough data points for Kolmogorv-Smirnov-test.')
elseif n == 4
    if nDnn >= 4
        p = 'sign. at alpha = 0.05';
    else
        p = 'not sign. at alpha = 0.05';
    end
elseif (n == 5) | (n == 5)
    if nDnn >= 5
        p = 'sign. at alpha = 0.05';
    else
        p = 'not sign. at alpha = 0.05';
    end
elseif (n >= 7) & (n <= 9)
    if nDnn >= 6
        p = 'sign. at alpha = 0.05';
    else
        p = 'not sign. at alpha = 0.05';
    end
elseif (n >= 10) & (n <= 13)
    if nDnn >= 7
        p = 'sign. at alpha = 0.05';
    else
        p = 'not sign. at alpha = 0.05';
    end
elseif (n >= 14) & (n <= 17)
    if nDnn >= 8
        p = 'sign. at alpha = 0.05';
    else
        p = 'not sign. at alpha = 0.05';
    end
elseif (n >= 18) & (n <= 22)
    if nDnn >= 9
        p = 'sign. at alpha = 0.05';
    else
        p = 'not sign. at alpha = 0.05';
    end
elseif (n >= 23) & (n <= 27)
    if nDnn >= 10
        p = 'sign. at alpha = 0.05';
    else
        p = 'not sign. at alpha = 0.05';
    end
elseif (n >= 28) & (n <= 32)
    if nDnn >= 11
        p = 'sign. at alpha = 0.05';
    else
        p = 'not sign. at alpha = 0.05';
    end
elseif (n >= 33) & (n <= 38)
    if nDnn >= 12
        p = 'sign. at alpha = 0.05';
    else
        p = 'not sign. at alpha = 0.05';
    end
else
    error('Data size out of range.')
end