function [ccr lags] = eegxcorr(eeg1,eeg2)
%EEGXCORR   Crosscorrelation fo EEG.
%   EEGXCORR(EEG1,EEG2) calculates crosscorrelogram for EEG
%   EEG1 and EEG2, using a +-50 ms time window.
%
%   [XC LAGS] = EEGXCORR(EEG1,EEG2) returns crosscorrelation in XC and a
%   corresponding vector of time lags in LAGS.
%
%   See also XCORR and CZXCORR.

% Input argument check
error(nargchk(2,2,nargin))

% Crosscorrelogram
sr = 1000;      % sampling rate
ccr = xcorr(eeg2,eeg1,0.05*sr);     % 1->2; window: -50 ms - 50 ms
lags = linspace(-50,50,length(ccr));

% Plot
if nargout == 0
    H1 = figure;
    bar(lags,ccr)
    set(gca,'XLim',[-50 50])
end