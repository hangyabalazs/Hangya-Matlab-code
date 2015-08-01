function czacorr(ncc)
%CZACORR   Autocorrelation.
%   CZACORR(VD) calculates autocorrelogram for discriminated unit VD, using
%   a +-200 ms time window.
%
%   See also XCORR and CZXCORR.

% Calculate spike times in milliseconds
sr = 1000;
nc = ncc * sr;
nc(nc<0.5) = [];

% Autocorrelogram
zunit1 = zeros(1,length(round(nc))+5);
zunit1(round(nc)) = 1;
acr = xcorr(zunit1,0.2*sr);
acr(length(acr)/2+0.5) = [];
acr = reshape(acr,length(acr)/100,100);     % window: -200 ms - 200 ms
sacr = sum(acr);

% Plot result
figure;
bar(linspace(-200,200,length(sacr)),sacr)