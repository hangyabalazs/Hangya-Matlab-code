function [T,seglen] = b_thres4
%THRES4   Automatic unit thersholding.
%   THRES4 determines the threshold through the follwing algorithm: it seeks the
%   upper limit of the noise zone and takes the minimum of the 16th to 20th highest
%   amplitudes (in order to exclude arteficial high values). Threshold is at the one
%   fourth distance of the above two values.
%
%   THRES4, unlike THRES segments long data and returns an array of threshold values.
%   [T,SEGLEN] = THRES4 returns the length of unit segments as well.
%
%   THRES4 can handle pseudounit as well. It does not plot the result.
%
%   See also THRES, THRES2, THRES3 and THRESRUN.

% Input arguments check
error(nargchk(0,0,nargin));

% Reading the input data
global IN
if isempty(IN)
    error('Data has not been inported.')
end
b_var2ws('in','caller')

% Assign segment length
seglen = 100000;

% Threshold pseudounit
if length(find(unit==max(unit))) > 3
    T = [];
    ind2 = 0;
    lenu = length(unit);
    while ind2 < lenu
        ind2 = ind2 + seglen;
        T(end+1) = 0.25;
    end
    return
end

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
% H = figure;
% plot(time,unit)

% Segmenting the unit
MS = [];
MI = [];
T = [];
ind2 = 0;
lenu = length(unit);
while ind2 < lenu
    ind1 = ind2 + 1;
    ind2 = ind2 + seglen;
    unitsegment = unit(ind1:min(ind2,lenu));
    
% Action potential amplitudes
    unit2 = unitsegment;
    for k = 1:20
        m(k) = max(unit2);
        fnd = find(unit2==m(k));
        unit2(fnd(1)) = [];
    end
    MS(end+1) = min(m(16:20));
    
% Noise
    z = mean(unitsegment);
    m = max(unitsegment);
    npc = 1;
    nsup = 0.99;
    ninf = 0.98;
    next = 0;
    while (npc > nsup | npc < ninf) & next <= 200
        if npc > ninf
            m = m - (m - z) / 2;
        else
            m = m + (m - z) / 2;
        end
        f = find(unitsegment<m);
        npc = length(f) / length(unitsegment);
        next = next + 1;
        if next == 201
            str = [fname(1:6) ' : Maximum number of iterations reached in noise limit search.'];
            disp(str)
        end
    end
    MI(end+1) = m;
    
% Threshold
    T(end+1) = MI(end) + (MS(end) - MI(end)) / 4;

% Plot threshold
%     line([time(ind1) time(min(ind2,lenu))],[MI(end) MI(end)],'Color','r')
%     line([time(ind1) time(min(ind2,lenu))],[MS(end) MS(end)],'Color','r')
%     line([time(ind1) time(min(ind2,lenu))],[T(end) T(end)],'Color','g')
end
% t = ['THRESHOLD ',fname(1:3),' ',fname(5:6),' ',num2str(datinx1),' ',num2str(datinx2)];
% title(t);