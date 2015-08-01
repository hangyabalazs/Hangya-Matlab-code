function wave = b_wavethres(wave)
%WAVETHRES    Wavelet thresholder.
%   WAVETHRES(WAVE) finds the maximum of wavelet power WAVE and replaces the
%   values lower than 5 per cent of the maximum with zeros.
%
%   See also WAVELET.

% Input argument check
error(nargchk(1,1,nargin))

% Thresholding
M = max(max(wave));
CV = M / 20;    % cut value: 5% of the maximum
wave(find(wave<CV)) = 0;