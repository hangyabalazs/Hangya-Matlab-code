function [H1 H2 trsc] = somccg(ncc1,ncc2,wn)
%SOMCCG   Crosscorrelation.
%   SOMCCG(VD1,VD2) calculates crosscorrelogram for discriminated units
%   VD1 and VD2, using a +-50 ms time window. Crosscorrelogram is
%   normalized with a shuffled ISI crosscorrelogram to remove random
%   correlations. A significance level of p=0.0013 is indicated on the
%   result. Original and normalized crosscorrelograms are plotted.
%
%   [H1 H2] = SOMCCG(VD1,VD2) returns the handles of the resulting plots.
%
%   [H1 H2 TRSC] = SOMCCG(VD1,VD2) returns transmission success rate for
%   excitatory connections as well. Note, that in case one of the neurons
%   is an interneuron, the interneuron has to be the first input argument
%   and pyramidal cell has to be the second input argument for transmission
%   success rate calculation.
%
%   [H1 H2 TRSC] = SOMCCG(VD1,VD2,WN) uses WN input argument as time window
%   for the cross-correlogram. WN should be given in milliseconds.
%   Normalized crosscorrelogram will not be effected.
%
%   See also XCORR and CZXCORR.

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
ccr = xcorr(zunit2,zunit1,wn2*sr);     % 1->2; window: -50 ms - 50 ms
% ccr = ccr / length(nc2);     % norm. with the no. of ref. events to get transmission prob.
if isequal(ncc1,ncc2)
    ccr(length(ccr)/2+0.5) = [];    % auto-correlation: drop middle bin
end
H1 = figure;
bar(linspace(-wn,wn,length(ccr)),ccr,'FaceColor','black')
set(gca,'XLim',[-wn wn])

% ISI shuffle
for k = 1:10
    rnd = rand(1) * 100 + 50;
    shf = round(rnd);
    pzunit = [zunit2(shf+1:end) zunit2(1:shf)];
    
    % Random crosscorrelogram
    pccr = xcorr(pzunit,zunit1,0.05*sr);
    pccr = pccr / length(nc2);
    str = ['p' num2str(k) '=pccr;'];
    eval(str)
end
pccr = mean([p1;p2;p3;p4;p5;p6;p7;p8;p9;p10]);
% figure;
% bar(pccr)

% Normalized crosscorrelogram
if isequal(ncc1,ncc2)
    H2 = [];    % auto-correlation: return
    trsc = [];
    return
end
if ~isequal(wn,50)
    ccr = xcorr(zunit2,zunit1,0.05*sr);     % 1->2; window: -50 ms - 50 ms
    ccr = ccr / length(nc2);     % norm. with the no. of ref. events to get transmission prob.
end
nccr = ccr - pccr;
trsc = max(nccr(45:55));      % transmission success rate
H2 = figure;
bar(linspace(-50,50,length(nccr)),nccr,'FaceColor','black')
set(gca,'XLim',[-50 50])

mnc = mean(nccr);
sdc = std(nccr);
thr = mnc + 3 * sdc;
nthr = mnc - 3 * sdc;
line([-50 50],[thr thr],'Color','red')      % significance level: p=0.0013
line([-50 50],[nthr nthr],'Color','red')
line([-50 50],[nthr nthr],'Color','red')