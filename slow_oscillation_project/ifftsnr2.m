function ifftsnr2(pat,patno,eg,nm_rows,nm_cols,esg)
%IFFTSNR2   Fourier-based SNR calculation.
%   IFFTSNR2 calculates Fourier spectra and estimates signal-to-noise ratio
%   as the ratio of spectral power in the delta band relative to the 4-40
%   Hz frequency band. Channel amplitudes are also calculated and saved.
%
%   Note: prefiltering data (40 Hz, lowpass) is recommended (see IFILTER)!
%
%   IFFTSNR2 can be called by ICALLER2. Optional input arguments:
%       PAT: patient code and name
%       PATNO: patient code
%       EG: EEG file no.
%       NM_ROWS: number of grid rows
%       NM_COLS: number of grid columns
%
%   See also IFFTSNR, IMISHIFT and IFILTER3.

% Input argument check
error(nargchk(0,6,nargin))
if nargin <6
    esg = [1040 1100];
end
if nargin < 5
    nm_cols = 5;           % number of columns on the electrode grid
end
if nargin < 4
    nm_rows = 4;           % number of rows on the electrode grid
end
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
global DATAPATH
resdir = [DATAPATH 'Ulbert\OITI_' patno '_EEG_' eg '\SNR\'];
try
    cd(resdir)
catch
    mkdir(resdir)
    cd(resdir)
end
mm = pwd;
cd(resdir)

% Load raw data
global DATADIR
ddir = [DATADIR '\human_SO\' pat '\grid\mat_EEG_' eg '\'];
fnm = ['EEG_' eg '_' num2str(esg(1)) '_' num2str(esg(2)) '_rs_filt01_40.mat'];
load([ddir fnm])
sr = 1000;

% FFT
chno = size(data,2);    % prefiltering data (0.1-40 Hz) is recommended (see IFILTER3)!
nm_rows = 4;
if chno == 16
    nm_cols = 4;
else
    nm_cols = 5;
end
nfft = 8192;
EcogFFT = zeros(chno,nfft);
Amplitude = zeros(1,chno);
SNR = zeros(1,chno);
SNR2 = zeros(1,chno);
SNR3 = zeros(1,chno);
Noise = zeros(1,chno);
H = figure;
for x = 1:chno
    ch1 = data(:,x)';
    [EcogFFT(x,:),w] = b_fft2(ch1,sr,nfft);
    SNR(x) = mean(EcogFFT(x,(w>=0.5&w<=4))) / mean(EcogFFT(x,(w>4&w<=40)));    % proportion of delta power to non-delta power
    SNR2(x) = mean(EcogFFT(x,(w>=0.5&w<=4))) / mean(EcogFFT(x,(w>=0.5&w<=40)));   % proportion of delta power to all power
    SNR3(x) = max(EcogFFT(x,(w>=0.5&w<=4))) / mean(EcogFFT(x,(w>=0.5&w<=40)));   % max delta power to all power (Hurtado et al 2004)
    Noise(x) = mean(EcogFFT(x,(w>4&w<=40)));   % noise power
    ch1_norm = ch1 - mean(ch1);
    Amplitude(x) = mean(abs(ch1_norm));
    
    % Plot
    figure(H);
    subplot(nm_rows,nm_cols,x);
    plot(w(2:nfft/32),EcogFFT(x,2:nfft/32))
end

% Save
cd(resdir)
save(['SNR' eg '.mat'],'EcogFFT','SNR','SNR2','SNR3','Noise','Amplitude','w')
saveas(H,['SNR' eg '.fig'])
cd(mm)