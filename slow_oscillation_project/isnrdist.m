function [TIM SNR] = isnrdist(pat,patno,eg)
%SNRDIST  SNR distribution.
%   [TIM SNR] = ISNRDIST(PAT,PATNO,EG) returns numbers of propagating waves
%   (TIM) and signal-to-noise ratio (SNR) for all channels. Input
%   parameters:
%       PAT: patient ncode and name
%       PATNO: patient number
%       EG: EEG number
%
%   See also IFFTSNR and ISNRDIST_CALL.

% Input argument check
error(nargchk(0,3,nargin))
if nargin < 3
    eg = '46';
end
if nargin < 2
    patno = num2str('40');
end
if nargin < 1
    pat = 'oiti40_mojzsis';
end

% Directories
global DATADIR
global DATAPATH
bi = [DATADIR 'human_SO\' pat '\' pat '.jpg'];
inpdir = [DATAPATH 'Ulbert\OITI_' patno '_EEG_' eg '\MImap\'];
inpdir_snr = [DATAPATH 'Ulbert\OITI_' patno '_EEG_' eg '\SNR\'];

% Load
fn = [inpdir 'MIshiftfine.mat'];    % MI map
load(fn)
fn = [inpdir_snr 'SNR' eg '.mat'];    % MI map
load(fn)

% Significance level
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

% Typical propagation rates
TIM = nansum(logical(nansum(Adj,2)),3) / 10;
TIM = TIM';
TIM = TIM(:);

% Signal-to-noise ratio
SNR = Noise';
SNR = SNR(:);
if isequal(patno,'37')
    TIM(18) = [];
    SNR(18) = [];
end
if isequal(patno,'N2')
    TIM(5) = [];
    SNR(5) = [];
    TIM(17) = [];
    SNR(17) = [];
    TIM(18) = [];
    SNR(18) = [];
end
patno
max(SNR) / min(SNR)    