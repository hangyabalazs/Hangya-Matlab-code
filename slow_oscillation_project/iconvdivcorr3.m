function iconvdivcorr3(pno,pt,egs)
%ICONVDIVCORR3   Correlations of convergence and divergence patterns.
%   ICONVDIVCORR3 calculates correlation between convergence and divergence
%   pattern of the electrode receiving the most convergent input. See
%   IMISORUN for details on the above variables.
%
%   ICONVDIVCORR3(PNO,PT,EGS) requires 3 string input arguments:
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
wf = cell(1,legs);
wt = cell(1,legs);
for sgs = 1:legs
    eg = egs{sgs};
    inpdir = [DATAPATH 'Ulbert\OITI_' pno '_EEG_' eg '\Mtx\'];
    inpdir2 = [DATAPATH 'Ulbert\OITI_' pno '_EEG_' eg '\MImap\'];
    
    fn = [inpdir 'mtxs_' eg '.mat'];    % MI map
    load(fn)
    nm_rows = size(Convergence,1);
    nm_cols = size(Convergence,2);
    [Vfrom Vto Vend nm_arrows] = inptrf(inpdir2,eg,nm_rows,nm_cols);
    
    sCnv = sum(Convergence,3);
    if sgs == 1
        mxc = find(sCnv'==max(sCnv(:)));
    end
    chnum = nm_rows * nm_cols;
    wf{sgs} = zeros(1,chnum);   % incoming arrows to the electrode with the largest no. of converging arrows
    for k = 1:nm_arrows
        if Vto(k) == mxc && any(Vto==mxc&((Vend>Vend(k)&Vend<Vend(k)+100)|(Vend<Vend(k)&Vend>Vend(k)-100)))
            ltr = Vfrom(k);
            wf{sgs}(ltr) = wf{sgs}(ltr) + 1;
        end
    end
    
    sDiv = sum(Divergence,3);
    chnum = nm_rows * nm_cols;
    wt{sgs} = zeros(1,chnum);   % incoming arrows to the electrode with the largest no. of converging arrows
    for k = 1:nm_arrows
        if Vfrom(k) == mxc
            ltr = Vto(k);
            wt{sgs}(ltr) = wt{sgs}(ltr) + 1;
        end
    end
end

% Correlations
next = 1;
R = zeros(1,legs);
p = zeros(1,legs);
for xi = 1:legs
    x = wf{xi};
    y = wt{xi};
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
xlsname = [resdir 'maxconvdivcorr.xls'];   % write results to excel file
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
function [Vfrom Vto Vend nm_arrows] = inptrf(inpdir,eg,nm_rows,nm_cols)

% Original MI map
fn = [inpdir 'MIshiftfine.mat'];    % MI map
load(fn)

% Significance levels
fn = [inpdir 'siglev_EEG' eg];
load(fn)
siglev = sl(4,2);   % sig. lev.: 0.0001

% Transform input variables
Adj = rIMax;        % adjacency matrix
Adj(rIMax<siglev) = NaN;
Adj(isnan(Adj)) = 0;
Adj(Adj>0) = 1;     % time-varying adjacency matrix
rIM = rIMax;
rIM(rIMax<siglev) = NaN;
rIM(isnan(rIM)) = 0;
rIML = rIMaxLoc;
rIML(rIMax<siglev) = NaN;
rIML(isnan(rIML)) = 0;

% Recording properties
chnum = nm_rows * nm_cols;   % number of grid electrodes
tmax = size(Adj,3);    % time axis length

% Arrow starts and ends
Astart = Adj .* permute(repmat(1:100:tmax*100,[chnum 1 chnum]),[1 3 2]);
Aend = Astart + rIML;
[Vfrom Vto Vs] = ind2sub(size(Adj),find(Adj));     % vector indeces
Vstart = (Vs - 1) * 100 + 1;
Vend = Vstart + rIML(find(Adj));
[xV yV zV] = ind2sub(size(Adj),find(Adj));     % vector indeces
nm_arrows = length(xV);    % number of vectors