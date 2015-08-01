function snwphase2

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

% Respiration phase
% resp_phase = angle(hilbert(standardize(fsniff)));    % Hilbert-transform
% figure
% plot((sniff-mean(sniff))/std(sniff))
% hold on
% plot((fsniff-mean(fsniff))/std(fsniff),'r')
% plot(resp_phase,'c')
% figure
% plot((sniff-mean(sniff))/std(sniff))
% hold on
% plot(resp_phase,'r')

% Discriminate whisking
vdisc = b_udisc(fwhisk');

% Respiration phase - whisking relationship
[hang_sniff hmvl_sniff ftm_sniff bahee_sniff ahee_sniff inx_sniff H_sniff Hc_sniff] = ...
    phasehist(sniff,vdisc,sr,eps,30,1024,0,0);

% EEG phase - whisking relationship
[hang_eeg hmvl_eeg ftm_eeg bahee_eeg ahee_eeg inx_eeg H_eeg Hc_eeg] = ...
    phasehist(eeg,vdisc,sr,eps,30,1024,0,0);

% Phase of out-of-sniffing-phase whisking rel. to EEG
bh = bahee_sniff - hang_sniff;  % centralized sniffing phase
bh2 = bh;
bh2(bh<0) = bh2(bh<0) + 2 * pi;
bh2(bh>=0) = bh2(bh>=0) - 2 * pi;
opi = min(abs([bh'; bh2']));    % phase difference rel. to mean phase
op = opi < pi / 2;      % optimal phase indeces (logical)
oop = ~op;      % suboptimal phase indeces (logical)
vdisc2 = vdisc;
vdisc2(inx_sniff) = [];
[hang_eeg hmvl_eeg ftm_eeg bahee_eeg ahee_eeg H_eeg Hc_eeg] = ...
    phasehist(eeg,vdisc2(oop),sr,eps,30,1024,0,0);      % phases rel. to EEG of whisking with suboptimal phase rel. to sniffing
[hang_eeg hmvl_eeg ftm_eeg bahee_eeg ahee_eeg H_eeg Hc_eeg] = ...
    phasehist(sniff,vdisc2(oop),sr,eps,30,1024,0,0);    % control: should be the hist of non-optimal phases; inx should be empty
[hang_eeg hmvl_eeg ftm_eeg bahee_eeg ahee_eeg H_eeg Hc_eeg] = ...
    phasehist(eeg,vdisc2(op),sr,eps,30,1024,0,0);


% Root mean square of whisking
lenw = length(fwhisk);
wl = sr / 100;     % 10 ms windows (for 1 kHz samp. rate)
lle = floor(lenw/wl) * wl;
wh2 = reshape(fwhisk(1:lle),wl,lle/wl);
rms_whisk = sqrt(sum(wh2.^2)) / sqrt(wl);
rms_whisk = rms_whisk';
resp_phase2 = resp_phase(1:lle);
resp_phase2 = resp_phase2(round(wl/2):wl:end);   % downsampled resp. phase aligned to whisking RMS



