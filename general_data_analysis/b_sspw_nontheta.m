function [bg,fin] = b_sspw_nontheta(eeg)
%SSPW_NONTHETA   Sharpwave detection in non-theta segments.
%   SSPW_NONTHETA uses the following double criterium for sharpwaves: root-mean-square
%   of EEG under sharpwave should reach mean(RMS) + std(RMS) and peak RMS of sharpwave 
%   should reach mean(RMS) + 5 * std(RMS).
%
%   Outputs: first and last points of sharpwaves.
%
%   Note: SHORT_WAVEPHASE is commented out!
%
%   See also NONTHETA.

% Input argument check
error(nargchk(1,1,nargin))

% Filtering
bg = [];
fin = [];
b = fir1(2048,[0.018 0.028]);
hfeeg = filtfilt(b,1,eeg);
ln = fix(length(hfeeg)/100);

% Calculating RMS (root-mean-square)
rmsw = zeros(1,ln);
for i = 0:ln-1
    rmsw(i+1) = norm(hfeeg((1+i*100):(100+i*100))) / 100;
end
rmseeg = mean(rmsw);
sdrms = std(rmsw);
spwprop = rmsw / sdrms;

% Detection of sharpwaves (5SD criterium)
sselector = zeros(1,ln);
nss = length(rmsw);
for i = 1:nss
    if rmsw(1,i) > rmseeg + sdrms   % first criterium: RMS under sharpwave should reach
        sselector(1,i) = 1;         % mean(RMS) + std(RMS)
    end
end
spws = hfeeg;
dss = diff(sselector);
dss1 = zeros(1,nss-1);
for i = 1:nss-1
    if dss(1,i) ~= 0
        dss1(1,i) = 1;
    end
end
if sselector(1) == 1
    dss1(1) = 1;
end
if sselector(end) == 1
    dss1(end) = 1;
end
dss2 = find(dss1(1,:)==1);
dss3 = dss2;
dss4 = dss3 .* 100;
dss5 = dss4(1,1:2:end);
dss6 = dss4(1,2:2:end);
if length(dss5) > length(dss6)
    dss5 = dss5(1:length(dss5)-1);
end
if length(dss6) > length(dss5)
    dss6 = dss6(1:length(dss6)-1);
end
dss7 = dss6 - dss5;
peaks = zeros(1,length(dss7));      % finding sharpwave peaks
for i = 1:length(dss7)
    peaks(i) = find(rmsw(dss5(i)/100:dss6(i)/100)==max(rmsw(dss5(i)/100:dss6(i)/100)));
end
peaks = peaks - 1;
fpeaks = peaks * 100;
fpeaks1 = dss5 + fpeaks;
trmsw = rmsw((dss5/100)+peaks);
speaks = zeros(1,length(trmsw));
for i = 1:length(trmsw)
    if trmsw(i) > rmseeg + 5*sdrms  % second criterium: peak should reach mean(RMS) + 5 * std(RMS)
        speaks(i) = 1;
    end
end
pos = find(speaks~=0);
spwn = length(pos);
bg = dss5(pos);     % first points of sharp waves in data points
fin = dss6(pos);    % last points of sharp waves in data points