function snwcoh

load('C:\Balazs\_data\HSW\R3\2008-07-09_14-46-47_R3\ALIGNED_SIGNALS.mat')

% Linearize data matrices 
PELLET_H1 = -PELLET_H1';
PELLET_SN = PELLET_SN';
PELLET_W1 = -PELLET_W1';
eeg = PELLET_H1(:);
sniff = PELLET_SN(:);
whisk = PELLET_W1(:);

% Resample
sr = 1 / (ALIGNED_TIME(2) - ALIGNED_TIME(1));
[p q] = rat(1000/sr);
eeg = resample(eeg,p,q);    % resample at 1000 Hz
sniff = resample(sniff,p,q);
whisk = resample(whisk,p,q);
sr = 1000;

% Filter
nqf = sr / 2;
% flt = fir1(1024,30/nqf,'low');
% fsniff = filtfilt(flt,1,sniff);     % lowpass filter sniffing at 30 Hz
flt = fir1(4096,10/nqf,'high');
fwhisk = filtfilt(flt,1,whisk);     % highpass filter whisking at 10 Hz
% flt = fir1(1024,[4 12]/nqf,'band');     % filter EEG in the theta band
% feeg = filtfilt(flt,1,eeg);
% figure
% plot(whisk)
% hold on
% plot(fwhisk,'r')
% figure
% plot(eeg)
% hold on
% plot(feeg,'r')

% Root mean square of whisking
lenw = length(fwhisk);
wl = sr / 100;     % 10 ms windows (for 1 kHz samp. rate)
lle = floor(lenw/wl) * wl;
wh2 = reshape(fwhisk(1:lle),wl,lle/wl);
rms_whisk = sqrt(sum(wh2.^2)) / sqrt(wl);
rms_whisk = rms_whisk';
sniff2 = sniff(1:lle);
sniff2 = sniff2(round(wl/2):wl:end);   % downsampled sniff aligned to whisking RMS
eeg2 = eeg(1:lle);
eeg2 = eeg2(round(wl/2):wl:end);   % downsampled EEG aligned to whisking RMS

% Coherence
figure
mscohere(fwhisk,eeg,hanning(1024),512,1024,sr)
figure
mscohere(fwhisk,sniff,hanning(1024),512,1024,sr)
figure
mscohere(rms_whisk,eeg2,hanning(1024),512,1024,sr/wl)
figure
mscohere(rms_whisk,sniff2,hanning(1024),512,1024,sr/wl)
1