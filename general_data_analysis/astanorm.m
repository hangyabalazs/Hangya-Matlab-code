function [sta sta_index1 sta_index2 lenv] = astanorm(vdisc,eeg,wn)
%ASTANORM    Normalized Spike Triggered Average.
%   [S O1 O2 L] = ASTANORM(VD,EEG,WN) calculates STA of EEG and discriminated
%   unit (VD) using WN windowsize. Each EEG window is standardized before
%   average calculation. STA, maximum STA and maximum STA minus mean STA is 
%   returned in S, O1 and O2. Number of spikes is returned in L.

% Standardize EEG
eeg = (eeg - mean(eeg)) / std(eeg);

% Calculate STA
wn2 = round(wn/2);
vdisc = vdisc(find(vdisc-wn2>0&vdisc+wn2<=length(eeg)));
lenv = length(vdisc);
st = zeros(lenv,2*wn2+1);
for t = 1:lenv
    eeg2 = eeg(vdisc(t)-wn2:vdisc(t)+wn2);
    st(t,1:2*wn2+1) = eeg2;
end
sta = mean(st,1);

% Output
sta_index1 = max(sta) - mean(sta);
sta_index2 = max(sta);