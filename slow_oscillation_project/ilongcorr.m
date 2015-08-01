function ilongcorr
%ILONGCORR   Correlations for hubs.
%   ILONGCORR calculates convergence and divergence correlations for long
%   and short segments. Results are saved in an Excel file.
%
%   See also IWITHINSUBJCORR.

% Directories
global DATAPATH
pno = '31';
pt = 'virag';
eg{1} = '285';
eg{2} = '285_long';
pat = ['oiti' pno '_' pt];
resdir = [DATAPATH 'Ulbert\Summary_resubmission\'];

% Load
sCnv = cell(1,2);
sDiv = cell(1,2);
for sgs = 1:2
    inpdir = [DATAPATH 'Ulbert\OITI_' pno '_EEG_' eg{sgs} '\Mtx\'];
    
    fn = [inpdir 'mtxs_' eg{sgs} '.mat'];    % MI map
    load(fn)
    
    sCnv{sgs} = sum(Convergence,3);
    sDiv{sgs} = sum(Divergence,3);
end

% Correlations - convergence
x = sCnv{1};
y = sCnv{2};
x = x(:);
y = y(:);
[b,bint,r,rint,stats] = regress(y,[ones(length(x),1) x]);
R1 = sqrt(stats(1));         % correlation coefficient (R-value of the regression)
preR = corrcoef(x,y);
R = preR(2);
F = stats(2);           % F-test for H0: all coeff.-s are zero
p = stats(3);           % F-test significance
xlsname = [resdir 'long_convergence.xls'];   % write results to excel file
if b_isfilename(xlsname)
    [ntz mtz atz] = xlsread(xlsname,'sheet1');
    pref = size(atz,1) + 1;
else
    pref = 1;
end
xlswrite(xlsname,{pat},'sheet1',['A' num2str(pref)])
xlswrite(xlsname,R','sheet1',['B' num2str(pref)])
xlswrite(xlsname,p','sheet1',['C' num2str(pref)])

% Correlations - divergence
x = sDiv{1};
y = sDiv{2};
x = x(:);
y = y(:);
[b,bint,r,rint,stats] = regress(y,[ones(length(x),1) x]);
R1 = sqrt(stats(1));         % correlation coefficient (R-value of the regression)
preR = corrcoef(x,y);
R = preR(2);
F = stats(2);           % F-test for H0: all coeff.-s are zero
p = stats(3);           % F-test significance
xlsname = [resdir 'long_divergence.xls'];   % write results to excel file
if b_isfilename(xlsname)
    [ntz mtz atz] = xlsread(xlsname,'sheet1');
    pref = size(atz,1) + 1;
else
    pref = 1;
end
xlswrite(xlsname,{pat},'sheet1',['A' num2str(pref)])
xlswrite(xlsname,R','sheet1',['B' num2str(pref)])
xlswrite(xlsname,p','sheet1',['C' num2str(pref)])