function H = b_thetaselector_beta2(wave,f,newstep)
%THETASELECTOR_BETA2    Theta and delta power ratio.
%   H = THETASELECTOR_BETA2(WAVE,F,NEWSTEP) calculates theta and delta power
%   ratio of WAVE power using F scalevector and NEWSTEP constant (which 
%   characterizes downsampling of the EEG sampled on 10 kHz). It plots the 
%   result on figure H.
%
%   Theta band: 2.5 - 6 Hz. Delta band: 1 - 2.5 Hz.
%
%   See also THETASELECTOR_BETA and THETASELECTOR_BETA3.

% Input arguments check
error(nargchk(3,3,nargin));

% Open figure window
close all
H = figure;
hold on

% Theta-delta power ratio
power = (abs(wave)) .^ 2;
clear wave
sw2 = size(power,2);
time = (([1:sw2] - 1) * newstep + 1) / 10000;
time = round(time*100) / 100;
thetaperdelta(power,f,'r')
b_retimefig(time)

% --------------------------------------------------------------------------------
function thetaperdelta(power,f,c)

% Computing theta power
fnd = find(f>6);
pwind_theta1 = fnd(end);
fnd = find(f<2.5);
pwind_theta2 = fnd(1);
thetapower = power(pwind_theta1:pwind_theta2,:);
sumthetapower = sum(thetapower);
clear thetapower

% Computing delta power
fnd = find(f>2.5);
pwind_delta1 = fnd(end);
fnd = find(f<1);
pwind_delta2 = fnd(1);
deltapower = power(pwind_delta1:pwind_delta2,:);
sumdeltapower = sum(deltapower);
clear deltapower

% Plot ratio
warning off
ratio = sumthetapower ./ sumdeltapower;
warning backtrace
plot(ratio,'Color',c)
y_lim = ylim;
x_lim = xlim;
axis([x_lim(1) length(ratio) y_lim(1) y_lim(2)]);
tt = '\theta power / \delta power';
title(tt)

% ----------------------------------------------------------------------------------
function b_retimefig(time)      %changes time labels to seconds
a0 = get(gca,'XTick');
iszero = 0;
if a0(1) == 0
    a0(1) = 1;
    iszero = 1;
end
a1 = time(a0);
if iszero
    a1(1) = 0;
end
set(gca,'XTickLabel',num2str(a1'))