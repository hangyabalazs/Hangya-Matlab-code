function [T,seglen] = b_thres3
%THRES3   Automatic unit thersholding.
%   THRES3 determines the threshold through the follwing algorithm: it seeks the
%   upper limit of the noise zone and takes the minimum of the 6th to 10th highest
%   amplitudes (in order to exclude arteficial high values). Threshold is at the one
%   third distance of the above two values.
%
%   THRES3, unlike THRES segments long data and returns an array of threshold values.
%   [T,SEGLEN] = THRES3 returns the length of unit segments as well.
%
%   THRES3, unlike THRES2 plots the result.
%
%   See also THRES, THRES2, THRES4 and THRESRUN.

% Input arguments check
error(nargchk(0,0,nargin));

% Reading the input data
global IN
if isempty(IN)
    error('Data has not been inported.')
end
b_var2ws('in','caller')

% Checking if base line is stable
d = 10000;
dd = 1;
mn = [];
while dd < time(end)
    mn(end+1) = mean(unit(dd:dd+d));
    dd = dd + d;
end
mn(end+1) = mean(unit(dd:end));
if std(mn) > 0.1
    str = [fname(1:6) ' : Instable base line.'];
    disp(str)
end

% Plot unit
H = figure;
plot(time,unit)

% Segmenting the unit
MS = [];
MI = [];
T = [];
seglen = 200000;
ind2 = 0;
lenu = length(unit);
while ind2 < lenu
    ind1 = ind2 + 1;
    ind2 = ind2 + seglen;
    unitsegment = unit(ind1:min(ind2,lenu));
    
% Action potential amplitudes
    unit2 = unitsegment;
    for k = 1:10
        m(k) = max(unit2);
        unit2(find(unit2==m(k))) = [];
    end
    MS(end+1) = min(m(6:10));
    
% Noise
    z = mean(unitsegment);
    m = max(unitsegment);
    npc = 1;
    nsup = 0.99;
    ninf = 0.98;
    while npc > nsup | npc < ninf
        if npc > ninf
            m = m - (m - z) / 2;
        else
            m = m + (m - z) / 2;
        end
        f = find(unitsegment<m);
        npc = length(f) / length(unitsegment);
    end
    MI(end+1) = m;
    
% Threshold
    T(end+1) = MI(end) + (MS(end) - MI(end)) / 3;

% Plot threshold
    line([time(ind1) time(min(ind2,lenu))],[MI(end) MI(end)],'Color','r')
    line([time(ind1) time(min(ind2,lenu))],[MS(end) MS(end)],'Color','r')
    line([time(ind1) time(min(ind2,lenu))],[T(end) T(end)],'Color','g')
end
t = ['THRESHOLD ',fname(1:3),' ',fname(5:6),' ',num2str(datinx1),' ',num2str(datinx2)];
title(t);