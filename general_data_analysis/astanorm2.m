function [sta sta_cv sta_index1 sta_index2 sta_amp lenv st] = astanorm2(vdisc,eeg,sr,wn,isnorm,nc)
%ASTANORMMOD    Normalized Spike Triggered Average.
%   [S CV O1 O2 A L] = ASTANORM(VD,EEG,SR,WN,NC,ISNORM) calculates STA of EEG
%   and discriminated unit. Input arguments:
%       VD: discriminated unit (in data points)
%       EEG: contiuous signal (EEG, LFP, ECoG)
%       SR: sampling rate
%       WN (optional): window size (default: 4 sec)
%       ISNORM (optional): 1 if normalized, 0 if not (default: 1)
%       NC (optional): normalizing constant (default: std(EEG))
%    Output arguments:
%       STA: normalized spike triggered average
%       CV: STA coefficient of variation (mean / SD)
%       O1: maximum STA
%       O2: maximum STA minus mean STA
%       A: STA amplitude
%       L: number of spikes used for the calculation
%       ST: STA stack for STA image
%
%   See also ASTANORM, RJSTA3 and STAFIG.

% Input argument check
error(nargchk(3,6,nargin))
if nargin < 6
    nc = std(eeg);
end
if nargin < 5
    isnorm = 1;
end
if nargin < 4
    wn = 4 * sr;
end

% Standardize EEG
if isnorm
    eeg = (eeg - mean(eeg)) / nc;
end

% Calculate STA
wn2 = round(wn/2);
vdisc = vdisc(vdisc-wn2>0&vdisc+wn2<=length(eeg));
lenv = length(vdisc);
eeginx = repmat(vdisc(:)-wn2,1,2*wn2+1) + repmat(0:2*wn2,lenv,1);
st = eeg(eeginx);
% for t = 1:lenv    % previous, slow way of calculating
%     eeg2 = eeg(vdisc(t)-wn2:vdisc(t)+wn2);
%     st(t,1:2*wn2+1) = eeg2;
% end
sta = mean(st,1);
sta_cv = sta ./ std(st,[],1);

% Output
sta_index1 = max(sta) - mean(sta);
sta_index2 = max(sta);
sta_amp = max(sta) - min(sta);