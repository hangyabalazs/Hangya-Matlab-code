function [H1 ccr lwr upr rccg] = somccg_conf_filter(ncc1,ncc2,wn,shuff_min,shuff_max,cutoff)
%SOMCCG_CONF_FILTER   Cross-correlation.
%   SOMCCG_CONF_FILTER(VD1,VD2) calculates cross-correlogram for
%   discriminated units VD1 and VD2, using a +-50 ms time window.
%   Confidance intervals are calculated based on cross-correlations of
%   shifted data. SOMCCG_CONF_FILTER calculates confidence interval from
%   high-pass filtered CCG, adding the low-pass filtered part back
%   afterwards (see XCORR_WRAND_FILTER). Filtering is at CUTOFF Hz (see
%   below).
%
%   H = SOMCCG_CONF_FILTER(VD1,VD2) returns the handles of the resulting
%   plot.
%
%   H = SOMCCG_CONF_FILTER(VD1,VD2,WN) uses WN input argument as time
%   window for the cross-correlogram. WN should be given in milliseconds.
%   Normalized cross-correlogram will not be effected.
%
%   [H1 CCR LWR UPR RCCG] = SOMCCG_CONF_FILTER(VD1,VD2,WN) returns
%   cross-correlogram (CCR), lower and upper confidance intervals (LWR and
%   UPR) and shifted cross-correlations (RCCG).
%
%   [H1 CCR LWR UPR RCCG] = SOMCCG_CONF_FILTER(VD1,VD2,WN,SHUFF_MIN,SHUFF_MAX) 
%   accepts input arguments for minimal and maximal shifts for shuffled
%   crosscorrelations in ms (default: 100 ms and 5 s).
%
%   [H1 CCR LWR UPR RCCG] = SOMCCG_CONF_FILTER(VD1,VD2,WN,SHUFF_MIN,SHUFF_MAX,CUTOFF)
%   accepts a CUTOFF parameter for filtering CCGs (see XCORR_WRAND_FILTER;
%   default is 4 Hz).
%
%   See also SOMCCG, XCORR_WRAND_FILTER and XCORR.

% Input argument check
error(nargchk(2,6,nargin))
if nargin < 3
    wn = 50;    % window size in ms
end
if nargin < 4
    shuff_min = 100;     % 100 ms
end
if nargin < 5
    shuff_max = 5000;    % 5 s
end
if nargin < 6
    cutoff = 4;     % default filter cutoff: 4 Hz
end
sr = 1000;
shuff_min = shuff_min * sr / 1000;
shuff_max = shuff_max * sr / 1000;
num_shuff = 5000;

% Calculate spike times in milliseconds
nc1 = ncc1 * sr;
nc2 = ncc2 * sr;
mn = min(nc1(1),nc2(1));  % only relative spike times count; avoid out of memory
nc1 = nc1 - mn;
nc2 = nc2 - mn;
nc1(nc1<0.5) = [];  % drop spikes before 0.5 ms(!) to avoid error in line 39
nc2(nc2<0.5) = [];
wn2 = wn / 1000;    % window size in seconds

% Crosscorrelogram
zunit1 = zeros(1,round(max([nc1; nc2]))+5);
zunit2 = zunit1;
zunit1(round(nc1)) = 1;
zunit2(round(nc2)) = 1;
shf = round(rand(1,shuff_max)*num_shuff+shuff_min);
[ccr,lags,rccg] = xcorr_wrand_filter(sr,cutoff,zunit2,zunit1,wn2*sr,shf);     % 1->2; window: -wn ms - wn ms
% ccr = ccr / length(nc2);     % norm. with the no. of ref. events to get transmission prob.
if isequal(ncc1,ncc2)
    ccr(length(ccr)/2+0.5) = [];    % auto-correlation: drop middle bin
    rccg(:,end) = [];
end

% Plot
H1 = figure;
time = linspace(-wn,wn,length(ccr));
bar(time,ccr,'FaceColor','black')
set(gca,'XLim',[-wn wn])

% Plot confidence interval
hold on
[lc nc] = size(rccg);
ptc = ceil(lc*0.0005);
upr = zeros(1,nc);
lwr = zeros(1,nc);
for k = 1:nc
    sts = sort(rccg(:,k),'ascend');
    upr(k) = sts(end-ptc);
    lwr(k) = sts(ptc);
end
plot(time,upr,'Color',[0.7 0.7 0.7])
plot(time,lwr,'Color',[0.7 0.7 0.7])