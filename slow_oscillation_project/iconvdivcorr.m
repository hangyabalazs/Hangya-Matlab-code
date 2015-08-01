function iconvdivcorr(pno,pt,egs)
%ICONVDIVCORR   Correlations of convergence and divergence patterns.
%   ICONVDIVCORR calculates correlation of convergence and divergence
%   matrices. See IMISORUN for details on the above variables.
%
%   ICONVDIVCORR(PNO,PT,EGS) requires 3 string input arguments:
%       PNO: patient ID
%       PT: patient name
%       EGS: EEG file IDs
%
%   See also IMISORUN and IMISORUN_MTX.

% Directories
global DATADIR
global DATAPATH
% pno = num2str(39);
% % pno = 'n1';
% pt = 'gaal'
% egs = [{'10'} {'14'} {'40'}];
pat = ['oiti' pno '_' pt];
resdir = [DATAPATH 'Ulbert\Summary_resubmission\'];

% Load
legs = length(egs);
sCnv = cell(1,legs);
sDiv = cell(1,legs);
for sgs = 1:legs
    eg = egs{sgs};
    inpdir = [DATAPATH 'Ulbert\OITI_' pno '_EEG_' eg '\Mtx\'];
    
    fn = [inpdir 'mtxs_' eg '.mat'];    % MI map
    load(fn)
    
    sCnv{sgs} = sum(Convergence,3);
    sDiv{sgs} = sum(Divergence,3);
end

% Correlations - convergence
next = 1;
R = zeros(1,legs);
p = zeros(1,legs);
for xi = 1:legs
    x = sCnv{xi};
    y = sDiv{xi};
    x = x(:) / sum(x(:));
    y = y(:) / sum(y(:));
    [b,bint,r,rint,stats] = regress(y,[ones(length(x),1) x]);
    R1 = sqrt(stats(1));         % correlation coefficient (R-value of the regression)
    preR = corrcoef(x,y);
    R(next) = preR(2);
    F = stats(2);           % F-test for H0: all coeff.-s are zero
    p(next) = stats(3);           % F-test significance
    next = next + 1;
end
xlsname = [resdir 'convdivcorr.xls'];   % write results to excel file
if b_isfilename(xlsname)
    [ntz mtz atz] = xlsread(xlsname,'sheet1');
    pref = size(atz,1) + 1;
else
    pref = 1;
end
xlswrite(xlsname,{pat},'sheet1',['A' num2str(pref)])
xlswrite(xlsname,R','sheet1',['B' num2str(pref)])
xlswrite(xlsname,p','sheet1',['C' num2str(pref)])



% -------------------------------------------------------------------------
function [inx1 inx2] = gridind2sub(ind,nm_cols)

inx1 = floor((ind-1)/nm_cols) + 1;
inx2 = rem(ind,nm_cols);
if inx2 == 0
    inx2 = nm_cols;
end