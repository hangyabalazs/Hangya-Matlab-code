function [ang inx cyclen1 feeg ahee] = aphase_gamma(eeg,vdisc,sr)
%APHASE_GAMMA    Phase angles for unit relative to EEG.
%   [A I] = APHASE_GAMMA(EEG,VDISC,SR) calculates Hilbert phase angles (A)
%   for discriminated unit (VDISC) relative to standardized EEG filtered in
%   the beta band, when sampling frequency is given in SR. Cycles not
%   fulfilling the following 2 criteria are discarded: (i) EEG amp. higher
%   then 2SD; (ii) min. 50 ms length. Indices of disclosed spikes of vdisc
%   are returned in I.
%
%   See also HILBERT.

% Filtering EEG
nqf = sr / 2;
flt = fir1(4096,[7 20]/nqf,'band');      % bandpass filtering
feeg = filtfilt(flt,1,eeg);
feeg = (feeg - mean(feeg)) / std(feeg);

% Hilbert transformation
ahee = angle(hilbert(feeg));

% Check SWS criteria:
% 1. discard cicles with EEG amp. lower then 2SD
% 2. discard cicles shorter then 50 ms
fn = find(-diff(ahee)>2*pi-0.3);
cyclen1 = mean(diff(fn)) / sr * 1000;   % cycle length in ms
sd = std(feeg);
inx = find(vdisc<fn(1));
for k = 1:length(fn)-1
    seeg = feeg(fn(k):fn(k+1));
    axs = max(seeg) - min(seeg);
    sahee = ahee(fn(k):fn(k+1));
    if (axs < 2 * sd)  | (fn(k+1) - fn(k) < 0.050 * sr)
        inx = [inx find(vdisc>fn(k)&vdisc<fn(k+1))];
    end
end
inx = [inx find(vdisc>fn(end))];
vdisc(inx) = [];
ang = ahee(vdisc);