function H1 = somccg_conf(ncc1,ncc2,wn)
%SOMCCG_CONF   Crosscorrelation.
%   SOMCCG_CONF(VD1,VD2) calculates crosscorrelogram for discriminated
%   units VD1 and VD2, using a +-50 ms time window. Confidance intervals
%   are calculated based on crosscorrelations of shifted data.
%
%   H = SOMCCG(VD1,VD2) returns the handles of the resulting plot.
%
%   H = SOMCCG(VD1,VD2,WN) uses WN input argument as time window for the
%   cross-correlogram. WN should be given in milliseconds. Normalized
%   crosscorrelogram will not be effected.
%
%   See also SOMCCG, XCORR and CZXCORR.

% Input argument check
error(nargchk(2,3,nargin))
if nargin < 3
    wn = 50;    % window size in ms
end

% Calculate spike times in milliseconds
sr = 1000;
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
shf = round(rand(1,1000)*1000+50);
[ccr,lags,rccg] = xcorr_wrand(zunit2,zunit1,wn2*sr,shf);     % 1->2; window: -200 ms - 200 ms
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
ptc = lc * 0.001 + 1;
upr = zeros(1,nc);
lwr = zeros(1,nc);
for k = 1:nc
    sts = sort(rccg(:,k),'ascend');
    upr(k) = sts(end-ptc);
    lwr(k) = sts(ptc);
end
plot(time,upr,'Color',[0.7 0.7 0.7])
plot(time,lwr,'Color',[0.7 0.7 0.7])