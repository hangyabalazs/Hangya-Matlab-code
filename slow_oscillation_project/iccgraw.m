function [rIxy,cIxy] = iccgraw(eeg1,eeg2,sr,overlap)
%ICCGRAW   Runs cross-correlation on a sequence of recordings.
%   ICCGRAW calculates cross-correlation (CCG) on 1000-point raw EEG
%   traces. It returns CCG for real data and shuffled segment controls.
%   Note, that input EEG should be downsampled on 1000 Hz!
%
%   ICCGRAW(EEG1,EEG2,SR,OLP) calculates CCG between EEG1 and EEG2 sampled
%   on SR using overlapping windows. Overlap is controlled by OLP (1 for
%   non-overlapping, 2 for 50% overlapping, etc.).
%
%   See also ICCG.

% Input argument check
error(nargchk(4,4,nargin))

% Create random eeg (segment shuffle)
seglen = length(eeg2) * sr / 1000;
ic1 = fix(seglen/sr);
increm = round(seglen/ic1);     % 1-1.2 sec. long segments
ed = {};
for t = 1:ic1
    ind1 = (t - 1) * increm + 1;
    ind2 = ind1 + increm - 1;
    ed{end+1} = eeg2(ind1:ind2);
end
led = length(ed);
rp = randperm(led);
while any(rp==[1:led])
    rp = randperm(led);
end
psed = [];
for j = 1:led
    psed = [psed ed{rp(j)}];
end

% Filter
eeg1 = lfilter(eeg1,1000);
eeg2 = lfilter(eeg2,1000);
psed = lfilter(psed,1000);

% Cross-correlation
elen = min(length(psed),seglen);
rIxy = ccg_line(eeg1(1:elen),eeg2(1:elen),1000,sr,overlap);    % real
cIxy = ccg_line(eeg1(1:elen),psed(1:elen),1000,sr,overlap);    % control



% -------------------------------------------------------------------------
% Cross-correlation
% -------------------------------------------------------------------------
function aIxy = ccg_line(data1,data2,WindowSize,sr,Overlap)

[k1 k2] = size(data1);

winlen = WindowSize;   % window size
maxi = floor(k2/winlen);
aIxy = [];

% CCG calculation
ovlp = Overlap;
for i = 1:maxi*ovlp-ovlp+1        % ABS LOOP
    inx1 = (i - 1) * winlen / ovlp + 1;  % Note: overlaping windows!
    inx1 = round(inx1);
    inx2 = inx1 + winlen - 1;

    y1 = data1(inx1:inx2);   % eeg1 segment
    y2 = data2(inx1:inx2);   % eeg2 segment
    
    Ixy = xcorr(y1,y2,sr);      % cross-correlation
    Ixy = max(Ixy(sr+1:end));
    
    aIxy = [aIxy Ixy];
end     % end of abs loop



% -------------------------------------------------------------------------
function data = lfilter(data,sr)
%IFILTER3   Filters multichannel EEG.
%   FD = IFILTER3(DATA,SR) filters DATA sampled on SR between 0.5 and 4 Hz
%   using a 4096 order lowpass FIR filter. Filtered data is returned in FD.
%
%   See also FIR1 and FILTFILT.

% Construct filter
nqf = sr / 2;      % Nyquist frequency
flt = fir1(2*4096,[0.5 4]/nqf,'band');      % bandpass filtering between 0.1 and 40 Hz

% Filter EEG
feeg = filtfilt(flt,1,data);

% Return filtered data
data = feeg;