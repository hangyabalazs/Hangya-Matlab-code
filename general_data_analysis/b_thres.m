function [T,H] = b_thres
%THRES  Automatic unit thersholding.
%   THRES determines the threshold through the follwing algorithm: it seeks the
%   upper limit of the noise zone and takes the minimum of the 6th to 10th highest
%   amplitudes (in order to exclude arteficial high values). Threshold is at the one
%   third distance of the above two values.
%
%   [T,H] = THRES returns the threshold T and handle H of the thresholded unit figure.
%
%   See also THRES2, THRES3, THRES4 and THRESRUN.

%Input arguments check
error(nargchk(0,0,nargin));

%Reading the input data
global IN
if isempty(IN)
    error('Data has not been inported.')
end
b_var2ws('in','caller')

% Checking if base line is stable
d = 10000;
dd = 1;
mn = [];
while dd + d < length(unit)
    mn(end+1) = mean(unit(dd:dd+d));
    dd = dd + d;
end
mn(end+1) = mean(unit(dd:end));
if std(mn) > 0.1
    str = [fname(1:6) ' : Instable base line.'];
    disp(str)
end

% Action potential amplitudes
unit2 = unit;
for k = 1:10
    m(k) = max(unit2);
    unit2(find(unit2==m(k))) = [];
end
MS = min(m(6:10));

% Noise
z = mean(unit);
m = max(unit);
npc = 1;
nsup = 0.96;
ninf = 0.95;
while npc > nsup | npc < ninf
    if npc > ninf
        m = m - (m - z) / 2;
    else
        m = m + (m - z) / 2;
    end
    f = find(unit<m);
    npc = length(f) / length(unit);
end
MI = m;

% Threshold
T = MI + (MS - MI) / 3;

% Plotting
H = figure;
plot(time,unit)
line([time(1) time(end)],[MI MI],'Color','r')
line([time(1) time(end)],[MS MS],'Color','r')
line([time(1) time(end)],[T T],'Color','g')
t = ['THRESHOLD ',fname(1:3),' ',fname(5:6),' ',num2str(datinx1),' ',num2str(datinx2)];
title(t);