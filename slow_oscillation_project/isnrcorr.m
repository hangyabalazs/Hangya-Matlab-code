function isnrcorr(pno,pt,eg)
%ISNRCORR   Correlations with SNR.
%   ISNRCORR calculates correlation between convergence/divergence and
%   different measures of signal-to-noise ratio (SNR) to test the
%   hypothesis whether differences in SNR contributes to differences in
%   network parameters.
%
%   ISNRCORR(PNO,PT,EG) requires 3 string input arguments:
%       PNO: patient ID
%       PT: patient name
%       EG: EEG file ID
%
%   See also ISNRCORR_CALL, IMISORUN, IMISORUN_MTX and IFFTSNR.

% Directories
global DATADIR
global DATAPATH
% pno = num2str(39);
% % pno = 'n1';
% pt = 'gaal'
% eg = '10';
pat = ['oiti' pno '_' pt];
inpdir = [DATAPATH 'Ulbert\OITI_' pno '_EEG_' eg '\SNR\'];
inpdir2 = [DATAPATH 'Ulbert\OITI_' pno '_EEG_' eg '\Mtx\'];
resdir = [DATAPATH 'Ulbert\SNR_Summary\'];

% Load
fn = [inpdir2 'mtxs_' eg '.mat'];    % MI map
load(fn)
fn = [inpdir 'SNR' eg '.mat'];    % MI map
load(fn)

% Correlations
SNR3
R = zeros(1,5);
p = zeros(1,5);
y = squeeze(sum(Convergence,3));
y = y';
y = y(:);   % to conform the order of channels set in imisorun
x = SNR';
[b,bint,r,rint,stats] = regress(y,[ones(length(x),1) x]);
R1 = sqrt(stats(1));         % correlation coefficient (R-value of the regression)
preR = corrcoef(x,y);
R(1) = preR(2);
F = stats(2);           % F-test for H0: all coeff.-s are zero
p(1) = stats(3);           % F-test significance

x = SNR2';
[b,bint,r,rint,stats] = regress(y,[ones(length(x),1) x]);
R1 = sqrt(stats(1));         % correlation coefficient (R-value of the regression)
preR = corrcoef(x,y);
R(2) = preR(2);
F = stats(2);           % F-test for H0: all coeff.-s are zero
p(2) = stats(3);           % F-test significance

x = SNR3';
[b,bint,r,rint,stats] = regress(y,[ones(length(x),1) x]);
R1 = sqrt(stats(1));         % correlation coefficient (R-value of the regression)
preR = corrcoef(x,y);
R(3) = preR(2);
F = stats(2);           % F-test for H0: all coeff.-s are zero
p(3) = stats(3);           % F-test significance

x = Noise';
[b,bint,r,rint,stats] = regress(y,[ones(length(x),1) x]);
R1 = sqrt(stats(1));         % correlation coefficient (R-value of the regression)
preR = corrcoef(x,y);
R(4) = preR(2);
F = stats(2);           % F-test for H0: all coeff.-s are zero
p(4) = stats(3);           % F-test significance

x = Amplitude';
[b,bint,r,rint,stats] = regress(y,[ones(length(x),1) x]);
R1 = sqrt(stats(1));         % correlation coefficient (R-value of the regression)
preR = corrcoef(x,y);
R(5) = preR(2);
F = stats(2);           % F-test for H0: all coeff.-s are zero
p(5) = stats(3);           % F-test significance

xlsname = [resdir 'SNRcorr_convergence2.xls'];   % write results to excel file
if b_isfilename(xlsname)
    [ntz mtz atz] = xlsread(xlsname,'sheet1');
    pref = size(atz,1) + 1;
else
    pref = 1;
end
% xlswrite(xlsname,{[pat ' ' eg]},'sheet1',['A' num2str(pref)])
% xlswrite(xlsname,R','sheet1',['B' num2str(pref)])
% xlswrite(xlsname,p','sheet1',['C' num2str(pref)])



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