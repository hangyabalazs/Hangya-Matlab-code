function [MaxLocIxy MaxIxy] = iccg(eeg1,eeg2,sr)
%ICCG   Linear cross-correlation.
%   [ML M] = ICCG(EEG1,EEG2,SR) calculates cross-correlation (CCG) on
%   1000-point non-overlapping raw EEG traces sampled on SR. Correlation is
%   calculated when shifting EEG2 by 0 - 1000 ms forward. Maximal corr. in
%   the function of the time shift is returned in M, maximum location is
%   returned in ML.
%   Note, that input EEG should be downsampled on 1000 Hz!
%
%   See also IMISHIFT.

% Input argument check
error(nargchk(3,3,nargin))

% MI significance levels based on shuffled segment EEG controls
siglev95 = 1.5958;      % 95% significance level
siglev99 = 1.9121;      % 99% significance level

% Cross-correlation calculation
leneeg = length(eeg1);
winlen = 1 * sr;   % window size: 1 sec.
maxi = floor(leneeg/winlen);
bno = fix(exp(0.626+0.4*log(winlen-1)));   % number of bins for distribution estimation
ovlp = 1;       % non-overlapping windows
maxii = maxi * ovlp - ovlp + 1;
miny1 = min(eeg1);
maxy1 = max(eeg1);
miny2 = min(eeg2);
maxy2 = max(eeg2);
MaxIxy = zeros(1,maxii);
MaxLocIxy = zeros(1,maxii);
for i = 1:maxii-1
    inx1 = (i - 1) * winlen / ovlp + 1;  % Note: overlaping windows if ovlp ~= 1
    inx1 = round(inx1);
    inx2 = inx1 + winlen - 1;
    
    Ixy = xcorr(eeg2(inx1:inx2),eeg1(inx1:inx2),sr);      % cross-correlation
    Ixy = Ixy(sr+1:end);
    MaxIxy(i) = max(Ixy);
    fnm = find(Ixy==MaxIxy(i));
    MaxLocIxy(i) = fnm(1);
    
end