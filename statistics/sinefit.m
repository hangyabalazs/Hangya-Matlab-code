function [x snss] = sinefit(data,sr,s,lb,ub)
%SINEFIT   Fit sine wave on data.
%   X = SINEFIT(DATA,SR,S,LP,UB) fits a sine wave on DATA sampled on SR
%   sampling rate using least square method. Start point (S), lower bound
%   (LB) and upper bound (UB) should contain three parameters: offset in
%   radians, frequency in Herz and amplitude. Output (X) is arranged in a
%   similar way.
%
%   [X SN] = SINEFIT(DATA,SR,S,LB,UB) returns the optimal sine wave in SN.
%
%   See also FMINSEARCHBND.

% Input argument check
[k1 k2] = size(data);
if isequal(k1,1)
    data = data';
end

% Define globals
global SINEFITDATA
SINEFITDATA = data;
global SINEFITSR
SINEFITSR = sr;

% Optimization
[x,f] = fminsearchbnd(@sinefiterf,s,lb,ub);
disp(['Estimation error: ' num2str(f)]);

% Plot result
sn = sin((0:(length(data)-1))/sr*2*pi*x(2)-x(1));
sns = sn / (max(sn) - min(sn)) * x(3);
snss = sns - mean(data);
figure
plot(data)
hold on
plot(snss,'r')

% Clear global
clear global SINEFITDATA
clear global SINEFITSR



% -------------------------------------------------------------------------
function fc = sinefiterf(xc)
%SINEFITERF   Error function for SINEFIT.
%   FC = FUNCT2(XC) returns least square error in FC. Input argument XC
%   should be 1-by-3 array containing offset, frequency and amplitude for
%   the sinewave to fit.
%
%   FMINSEARCHBND calls SINFITERF by minimization.
%
%   See also FMINSEARCHBND and SINEFIT.

% Get global variable
global SINEFITDATA
data = SINEFITDATA;
global SINEFITSR
sr = SINEFITSR;

% Define parameters
offset = xc(1);
fr = xc(2);
amp = xc(3);


% Calculate 'Rsquare'
n = length(data);
sn = sin((0:(n-1))/sr*2*pi*fr-offset);
sns = sn / (max(sn) - min(sn)) * amp;
snss = sns - mean(data);
Rsquare = sum((snss-data').^2);

% Output argument
fc = Rsquare;