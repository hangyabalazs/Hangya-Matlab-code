function iwithinsubjcorr(pno,pt,egs)
%IWITHINSUBJCORR   Within-subject correlations.
%   IWITHINSUBJCORR calculates within-subject correlation of convergence,
%   divergence, start point and first arrow matrices. See IMISORUN for
%   details on the above variables.
%
%   IWITHINSUBJCORR(PNO,PT,EGS) requires 3 string input arguments:
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
sStp = cell(1,legs);
sFia = cell(1,legs);
for sgs = 1:legs
    eg = egs{sgs};
    inpdir = [DATAPATH 'Ulbert\OITI_' pno '_EEG_' eg '\Mtx\'];
    
    fn = [inpdir 'mtxs_' eg '.mat'];    % MI map
    load(fn)
    
    sCnv{sgs} = sum(Convergence,3);
    sDiv{sgs} = sum(Divergence,3);
    sStp{sgs} = StartPoints;
    sFia{sgs} = FirstArrows;
end

% Correlations - convergence
next = 1;
R = zeros(1,legs);
p = zeros(1,legs);
for xi = 1:legs-1
    for yi = xi+1:legs
        x = sCnv{xi};
        y = sCnv{yi};
        x = x(:);
        y = y(:);
        [b,bint,r,rint,stats] = regress(y,[ones(length(x),1) x]);
        R1 = sqrt(stats(1));         % correlation coefficient (R-value of the regression)
        preR = corrcoef(x,y);
        R(next) = preR(2);
        F = stats(2);           % F-test for H0: all coeff.-s are zero
        p(next) = stats(3);           % F-test significance
        next = next + 1;
    end
end
xlsname = [resdir 'convergence.xls'];   % write results to excel file
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
next = 1;
R = zeros(1,legs);
p = zeros(1,legs);
for xi = 1:legs-1
    for yi = xi+1:legs
        x = sDiv{xi};
        y = sDiv{yi};
        x = x(:);
        y = y(:);
        [b,bint,r,rint,stats] = regress(y,[ones(length(x),1) x]);
        R1 = sqrt(stats(1));         % correlation coefficient (R-value of the regression)
        preR = corrcoef(x,y);
        R(next) = preR(2);
        F = stats(2);           % F-test for H0: all coeff.-s are zero
        p(next) = stats(3);           % F-test significance
        next = next + 1;
    end
end
xlsname = [resdir 'divergence.xls'];   % write results to excel file
if b_isfilename(xlsname)
    [ntz mtz atz] = xlsread(xlsname,'sheet1');
    pref = size(atz,1) + 1;
else
    pref = 1;
end
xlswrite(xlsname,{pat},'sheet1',['A' num2str(pref)])
xlswrite(xlsname,R','sheet1',['B' num2str(pref)])
xlswrite(xlsname,p','sheet1',['C' num2str(pref)])

% Correlations - startpoints
next = 1;
R = zeros(1,legs);
p = zeros(1,legs);
for xi = 1:legs-1
    for yi = xi+1:legs
        x = sStp{xi};
        y = sStp{yi};
        x = x(:);
        y = y(:);
        [b,bint,r,rint,stats] = regress(y,[ones(length(x),1) x]);
        R1 = sqrt(stats(1));         % correlation coefficient (R-value of the regression)
        preR = corrcoef(x,y);
        R(next) = preR(2);
        F = stats(2);           % F-test for H0: all coeff.-s are zero
        p(next) = stats(3);           % F-test significance
        next = next + 1;
    end
end
xlsname = [resdir 'startpoints.xls'];   % write results to excel file
if b_isfilename(xlsname)
    [ntz mtz atz] = xlsread(xlsname,'sheet1');
    pref = size(atz,1) + 1;
else
    pref = 1;
end
xlswrite(xlsname,{pat},'sheet1',['A' num2str(pref)])
xlswrite(xlsname,R','sheet1',['B' num2str(pref)])
xlswrite(xlsname,p','sheet1',['C' num2str(pref)])

% Correlations - firstarrows
next = 1;
R = zeros(1,legs);
p = zeros(1,legs);
for xi = 1:legs-1
    for yi = xi+1:legs
        x = sFia{xi};
        y = sFia{yi};
        x = x(:);
        y = y(:);
        [b,bint,r,rint,stats] = regress(y,[ones(length(x),1) x]);
        R1 = sqrt(stats(1));         % correlation coefficient (R-value of the regression)
        preR = corrcoef(x,y);
        R(next) = preR(2);
        F = stats(2);           % F-test for H0: all coeff.-s are zero
        p(next) = stats(3);           % F-test significance
        next = next + 1;
    end
end
xlsname = [resdir 'firstarrows.xls'];   % write results to excel file
if b_isfilename(xlsname)
    [ntz mtz atz] = xlsread(xlsname,'sheet1');
    pref = size(atz,1) + 1;
else
    pref = 1;
end
xlswrite(xlsname,{pat},'sheet1',['A' num2str(pref)])
xlswrite(xlsname,R','sheet1',['B' num2str(pref)])
xlswrite(xlsname,p','sheet1',['C' num2str(pref)])