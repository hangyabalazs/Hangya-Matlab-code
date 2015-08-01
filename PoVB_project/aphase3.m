function [ang cyclen1 cyclen2 cl phlen] = aphase3(eeg,vdisc,sr)
%APHASE3    Phase angles for unit relative to EEG.
%   [A,C1,C2,CL,PL] = APHASE2(EEG,VDISC,SR) calculates Hilbert phase angles
%   (A) for discriminated unit (VDISC) relative to EEG, when sampling 
%   frequency is given in SR. Cicles not fulfilling the following 2 criteria
%   are discarded: (i) EEG amp. higher then 2SD; (ii) min. 250 ms length. 
%   Mean cycle length before and after applying the criteria are returned 
%   in C1 and C2. CL cell contains length of the following cycles:
%       1. original
%       2. discarded upon 1st crit.
%       3. discarded upon 2nd crit.
%       4. all discarded
%       5. all remaining.
%   Cumulative 60 degree phase cycle segment length (in seconds) is
%   returned in PL (for phase firing rate calculation).
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
cl1 = [];
cl4 = [];
cl5 = [];
cl6 = [];
lens1 = 0;
lens2 = 0;
lens3 = 0;
lens4 = 0;
lens5 = 0;
lens6 = 0;
for k = 1:length(fn)-1
    seeg = feeg(fn(k):fn(k+1));
    axs = max(seeg) - min(seeg);
    sahee = ahee(fn(k):fn(k+1));
    if (axs < 2 * sd)  | (fn(k+1) - fn(k) < 0.25 * sr)
        inx = [inx find(vdisc>fn(k)&vdisc<fn(k+1))];
        cl5(end+1) = (fn(k+1) - fn(k)) / sr * 1000;   % discarded cycles' length in ms;
    else
        cl6(end+1) = (fn(k+1) - fn(k)) / sr * 1000;   % remaining cycles' length in ms;
    end
    if axs < 2 * sd
        cl1(end+1) = (fn(k+1) - fn(k)) / sr * 1000;   % cycle length in ms;
    end
    if fn(k+1) - fn(k) < 0.25 * sr
        cl4(end+1) = (fn(k+1) - fn(k)) / sr * 1000;   % cycle length in ms;
    end
    loclens1 = length(find(sahee>-pi&sahee<-pi+pi/3));
    lens1 = lens1 + loclens1;
    loclens2 = length(find(sahee>-pi+pi/3&sahee<-pi+2*pi/3));
    lens2 = lens2 + loclens2;
    loclens3 = length(find(sahee>-pi+2*pi/3&sahee<0));
    lens3 = lens3 + loclens3;
    loclens4 = length(find(sahee>0&sahee<pi/3));
    lens4 = lens4 + loclens4;
    loclens5 = length(find(sahee>pi/3&sahee<2*pi/3));
    lens5 = lens5 + loclens5;
    loclens6 = length(find(sahee>2*pi/3&sahee<pi));
    lens6 = lens6 + loclens6;
end
inx = [inx find(vdisc>fn(end))];
vdisc(inx) = [];
cyclen2 = mean(cl6) / sr * 1000;   % cycle length in ms
cl{1} = diff(fn) ./ sr .* 1000;   % original cycle length in ms
cl{2} = cl1;    % discarded upon 1st crit.
cl{3} = cl4;    % discarded upon 2nd crit.
cl{4} = cl5;    % all discarded
cl{5} = cl6;    % all remaining
ang = ahee(vdisc);
phlen = [lens1 lens2 lens3 lens4 lens5 lens6] / sr;     % cycle-segment length