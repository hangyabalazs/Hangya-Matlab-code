function [ang inx] = aphase4(eeg,vdisc,sr)
%APHASE4    Phase angles for unit relative to EEG.
%   [A I] = APHASE2(EEG,VDISC,SR) calculates Hilbert phase angles (A) for 
%   discriminated unit (VDISC) relative to EEG, when sampling frequency
%   is given in SR. Cycles not fulfilling the following 2 criteria are
%   discarded: (i) EEG amp. higher then 2SD; (ii) min. 250 ms length.
%   Indices of disclosed spikes of vdisc are returned in I.
%
%   See also HILBERT.

% Filtering EEG
nqf = sr / 2;
flt = fir1(4096,5/nqf,'low');      % lowpass filtering on 5 Hz
feeg = filtfilt(flt,1,eeg);

% Hilbert transformation
ahee = angle(hilbert(feeg));

% Check SWS criteria:
% 1. discard cicles with EEG amp. lower then 2SD
% 2. discard cicles shorter then 250 ms
fn = find(-diff(ahee)>2*pi-0.1);
cyclen1 = mean(diff(fn)) / sr * 1000;   % cycle length in ms
sd = std(feeg);
inx = find(vdisc<fn(1));
for k = 1:length(fn)-1
    seeg = feeg(fn(k):fn(k+1));
    axs = max(seeg) - min(seeg);
    sahee = ahee(fn(k):fn(k+1));
    if (axs < 2 * sd)  | (fn(k+1) - fn(k) < 0.25 * sr)
        inx = [inx find(vdisc>fn(k)&vdisc<fn(k+1))];
    end
end
inx = [inx find(vdisc>fn(end))];
vdisc(inx) = [];
ang = ahee(vdisc);