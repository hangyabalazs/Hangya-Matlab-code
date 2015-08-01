function data = ifilter3(data,sr)
%IFILTER3   Filters multichannel EEG.
%   FD = IFILTER3(DATA,SR) filters DATA sampled on SR between 0.1 and 40 Hz
%   using a 4096 order lowpass FIR filter. Filtered data is returned in FD.
%
%   See also FIR1 and FILTFILT.

% Construct filter
nqf = sr / 2;      % Nyquist frequency
flt = fir1(2*4096,[0.1 40]/nqf,'band');      % bandpass filtering between 0.1 and 40 Hz

% Filter EEG
chno = size(data,2);
feeg = zeros(size(data));
for x = 1:chno
    feeg(:,x) = filtfilt(flt,1,data(:,x));
end

% Return filtered data
data = feeg;