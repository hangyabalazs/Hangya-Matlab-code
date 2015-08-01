function data = ifilter(data,sr)
%IFILTER   Filters multichannel EEG.
%   FD = IFILTER(DATA,SR) filters DATA sampled on SR at 40 Hz using a 4096
%   order lowpass FIR filter. Filtered data is returned in FD.
%
%   See also FIR1 and FILTFILT.

% Construct filter
nqf = sr / 2;      % Nyquist frequency
flt = fir1(4096,40/nqf,'low');      % lowpass filtering at 40 Hz

% Filter EEG
chno = size(data,2);
feeg = zeros(size(data));
for x = 1:chno
    feeg(:,x) = filtfilt(flt,1,data(:,x));
end

% Return filtered data
data = feeg;