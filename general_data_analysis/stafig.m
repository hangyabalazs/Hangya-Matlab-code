function [H lwr upr] = stafig(sta,sta_rand,sta_amp,sta_amp_rand,sta_index2,sta_index2_rand,nn,wn,sr,titlestr)
%STAFIG   Plot function for Spike Triggered Average.
%   [H LWR UPR] = STAFIG(STA,STA_RAND,A,A_RAND,O2,O2_RAND,N,WN,SR,TITLESTR) plots STA.
%   Input arguments:
%       STA: Spike Triggered Average
%       STA_RAND: sample of randomized STA
%       A: SAT amplitude
%       A_RAND: amplitudes of random STA
%       O2: STA index2 = maximum STA minus mean STA
%       O2_RAND: STA index2 for random STA
%       N: number of spikes used for STA calculation
%       WN: window size
%       SR: sampling rate
%       TITLESTR: string input for figure title
%   Output parameter:
%       H: handle of the resulting figure
%       LWR: lower 95% confidence limit
%       UPR: upper 95% confidence limit
%   See also ASTANORM2 and RJSTA3.

% Plot
time = linspace(-wn/sr/2,wn/sr/2,length(sta));
H = figure;
plot(time,sta,'LineWidth',1.5)
hold on
ptc = length(sta_amp_rand) * 0.025 + 1;
upr = zeros(1,length(sta));
lwr = zeros(1,length(sta));

% Significance (2.5% on both sides)
for k = 1:size(sta,2)
    sts = sort(sta_rand(:,k),'ascend');
    upr(k) = sts(end-ptc);
    lwr(k) = sts(ptc);
end
plot(time,upr,'Color',[0.7 0.7 0.7])
plot(time,lwr,'Color',[0.7 0.7 0.7])

% Significance for STA amplitude
ach = allchild(H);
ax = findobj(ach,'type','axes');
title(ax(end),titlestr)
x_lim = xlim;
y_lim = ylim;
str = ['\it{Amp: }' '\bf ' num2str(sta_amp)];
if length(find(sta_amp>sta_amp_rand)) / length(sta_amp_rand) > 0.95
    clr = 'red';
else
    clr = 'black';
end
text(0.1,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color',clr)
str = ['\it{Max: }' '\bf ' num2str(sta_index2)];
if length(find(sta_index2>sta_index2_rand)) / length(sta_index2_rand) > 0.95
    clr = 'red';
else
    clr = 'black';
end
text(0.1,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color',clr')
str = ['\it{n: }' '\bf ' num2str(nn)];
text(0.1,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','black')