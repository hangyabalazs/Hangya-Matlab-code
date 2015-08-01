function H = b_thetaselector_beta(wave,f,newstep)
%THETASELECTOR_BETA    Theta and delta power ratio.
%   H = THETASELECTOR_BETA(WAVE,F,NEWSTEP) calculates theta and delta power
%   ratio of WAVE using F scalevector and NEWSTEP constant (which characterizes
%   downsampling of the EEG sampled on 10 kHz). It plots the result on figure H.
%   It also computes and plots theta-delta power ratio for contrasted wavelet 
%   power (so that low values do not distort the result).
%
%   Distribution of all power for normal and contrasted wavelet is plotted in 
%   another figure window.
%
%   Theta band: 2.5 - 6 Hz. Delta band: 1 - 2.5 Hz.
%
%   See also THETASELECTOR_BETA2 and THETASELECTOR_BETA3.

% Input arguments check
error(nargchk(3,3,nargin));

% Open figure window
close all
H = figure;
hold on

% Theta-delta power ratio
power = (abs(wave)) .^ 2;
clear wave
thetaperdelta(power,f,'r')

% Distribution of all power
sw1 = size(power,1);
sw2 = size(power,2);
time = (([1:sw2] - 1) * newstep + 1) / 10000;
lenw = sw1 * sw2;
linpower = reshape(power,1,lenw);
bno = 100;
B = figure;
[x y] = hist(linpower,bno);
bar(y,x);
title('Distribution of all power');

% Computations on contrasted wavelet
spower = sort(linpower);
clear linpower
lim1 = 0.2;
lim2 = 0.05;
step = -0.05;
legend_matrix{1} = ['percent: ' num2str('100')];
for percent = lim1:step:lim2
    lencicle = fix((lim2-lim1)/step+1);
    b = round((percent-lim1)/step+1);
    clr = (b - 1) / (lencicle - 1);
    c = [0 1-clr clr];
    [power,spower] = contrastwave(percent,power,spower,lenw,f,c,B,H);
    legend_matrix{b+1} = ['percent: ' num2str(percent*100)];
end
figure(H)
legend(legend_matrix)

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

% --------------------------------------------------------------------------------
function [power,spower] = contrastwave(percent,power,spower,lenw,f,c,B,H)
lenw_pcnt = round(percent * lenw);  % selecting the highest 15 per cent
spower = spower(end-lenw_pcnt+1:end);
reallimit = spower(1);
power(find(power<reallimit)) = 0;
figure(H)
thetaperdelta(power,f,c)
figure(B)
y_lim = ylim;
line([reallimit reallimit],[y_lim(1) y_lim(2)],'Color','red');
ycoord = y_lim(2) - (y_lim(2) - y_lim(1)) / 5;
text(reallimit,ycoord,num2str(percent));