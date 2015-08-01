function [st sta stats H] = stacall(vdisc,eeg,sr,wn)
%STACALL   Calls ASTANORM2 and STAFIG.
%   [ST STA STATS H] = STACALL(VDISC,EEG,SR,wn) calculates and plots Spike Triggered Average.
%   Input arguments:
%       VDISC: point process
%       EEG: continuous signal
%       SR: sampling rate
%       WN (optional): window size in data points (default: 4 sec)
%   Output argument:
%       ST: aligned data matrix
%       STA: spike triggered average
%       STATS: structure matrix of STA statistics with the followong fields:
%           CV: STA coefficient of variation (mean / SD)
%           STA_INDEX1: maximum STA
%           STA_INDEX2: maximum STA minus mean STA
%           STA_AMPLITUDE: STA amplitude
%           SPIKE_NUMBER: number of spikes used for the calculation
%           LOWER_LIMIT: lower 95% confidence limit
%           UPPER_LIMIT: upper 95% confidence interval
%       H: figure handle of the STA plot
%
%   See also ASTANORM2 and STAFIG.

% Input arguments for STA calculation
error(nargchk(3,4,nargin))
if nargin < 4
    wn = 4 * sr;
end
isnorm = 0;
wn2 = round(wn/2);
vdisc = vdisc(vdisc-wn2>0&vdisc+wn2<=length(eeg));
[sta sta_cv sta_index1 sta_index2 sta_amp lenv st] = astanorm2(vdisc,eeg,sr,wn,isnorm);
stats.sta_cv = sta_cv;
stats.sta_index1 = sta_index1;
stats.sta_index2 = sta_index2;
stats.sta_amplitude = sta_amp;
stats.spike_number = lenv;

% Randomized STA
nrs = 200;
sta_rand = zeros(nrs,wn+1);
wb = waitbar(0,'Generating surrogate data set. Please wait...','Name','Running STACALL...');  % progress indicator
global WB
WB(end+1) = wb;
for rpp = 1:nrs
    vdisc_rand = randpoisson(length(vdisc),length(eeg));
    [sta_rand(rpp,1:wn+1) temp_cv sta_index1_rand(rpp) sta_index2_rand(rpp)...
        sta_amp_rand(rpp)] = astanorm2(round(vdisc_rand),eeg,sr,wn,isnorm);
    waitbar(rpp/nrs)
end
close(wb)

% Plot
[H lwr upr] = stafig(sta,sta_rand,sta_amp,sta_amp_rand,sta_index2,sta_index2_rand,lenv,wn,sr,'');
stats.lower_limit = lwr;
stats.upper_limit = upr;
figure
imagesc(st)
colorbar