function [ang cyclen1 cyclen2 cl] = aphase(eeg,vdisc,sr)
%APHASE    Phase angles for unit relative to EEG.
%   [A,C1,C2,CL] = APHASE(EEG,VDISC,SR) calculates Hilbert phase angles (A)
%   for discriminated unit (VDISC) relative to EEG, when sampling frequency
%   is given in SR. Cicles not fulfilling the following 3 criteria are
%   discarded: (i) EEG amp. higher then 2SD; (ii) monotonic analytic
%   signal; (iii) min. 250 ms length. Mean cycle length before and after 
%   applying the criteria are returned in C1 and C2. CL cell contains 
%   length of the following cycles:
%       1. original
%       2. discarded upon 1st crit.
%       3. discarded upon 2nd crit.
%       4. discarded upon only 2nd crit.
%       5. discarded upon 3rd crit.
%       6. all discarded
%       7. all remaining.
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
% 2. discard cicles with non-monotonic analytic signal
% 3. discard cicles shorter then 250 ms
fn = find(-diff(ahee)>2*pi-0.1);
cyclen1 = mean(diff(fn)) / sr * 1000;   % cycle length in ms
sd = std(feeg);
inx = find(vdisc<fn(1));
cl1 = [];
cl2 = [];
cl3 = [];
cl4 = [];
cl5 = [];
cl6 = [];
for k = 1:length(fn)-1
    seeg = feeg(fn(k):fn(k+1));
    axs = max(seeg) - min(seeg);
    sahee = ahee(fn(k):fn(k+1));
    if (axs < 2 * sd)  | any(diff(sahee(2:end))<0) | (fn(k+1) - fn(k) < 0.25 * sr)
        inx = [inx find(vdisc>fn(k)&vdisc<fn(k+1))];
        cl5(end+1) = (fn(k+1) - fn(k)) / sr * 1000;   % discarded cycles' length in ms;
    else
        cl6(end+1) = (fn(k+1) - fn(k)) / sr * 1000;   % remaining cycles' length in ms;
    end
    if axs < 2 * sd
        cl1(end+1) = (fn(k+1) - fn(k)) / sr * 1000;   % cycle length in ms;
    end
    if any(diff(sahee(2:end))<0)
        cl2(end+1) = (fn(k+1) - fn(k)) / sr * 1000;   % cycle length in ms;
    end
    if any(diff(sahee(2:end))<0) & ~(axs < 2 * sd) & ~(fn(k+1) - fn(k) < 0.25 * sr)
        cl3(end+1) = (fn(k+1) - fn(k)) / sr * 1000;   % cycle length in ms;
    end
    if fn(k+1) - fn(k) < 0.25 * sr
        cl4(end+1) = (fn(k+1) - fn(k)) / sr * 1000;   % cycle length in ms;
    end
end
inx = [inx find(vdisc>fn(end))];
vdisc(inx) = [];
cyclen2 = mean(cl6) / sr * 1000;   % cycle length in ms
cl{1} = diff(fn) ./ sr .* 1000;   % original cycle length in ms
cl{2} = cl1;    % discarded upon 1st crit.
cl{3} = cl2;    % discarded upon 2nd crit.
cl{4} = cl3;    % discarded upon only 2nd crit.
cl{5} = cl4;    % discarded upon 3rd crit.
cl{6} = cl5;    % all discarded
cl{7} = cl6;    % all remaining
ang = ahee(vdisc);