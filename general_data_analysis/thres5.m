function [T,seglen] = thres5(unit,sr)
%THRES5   Automatic unit thersholding.
%   THRES5(UNIT,SR) determines the threshold for UNIT (sampled at SR Hz)
%   discrimination through the follwing algorithm: it seeks the upper limit
%   of the noise zone and takes the minimum of the 16th to 20th highest
%   amplitudes (in order to exclude arteficial high values). Threshold is
%   at the one third distance of the above two values.
%
%   THRES5, unlike THRES segments long data and returns an array of
%   threshold values. [T,SEGLEN] = THRES5(UNIT,SR) returns the length of
%   unit segments as well.
%
%   THRES5 can handle pseudounit as well. It does not plot the result.
%
%   See also THRES4 and THRESRUN2.

% Input arguments check
error(nargchk(2,2,nargin));

% Assign segment length
seglen = 20 * sr;   % 20 s segments

% Threshold pseudounit
lenu = length(unit);
if length(find(unit==max(unit))) > lenu / 2
    T = [];
    ind2 = 0;
    while ind2 < lenu
        ind2 = ind2 + seglen;
        T(end+1) = 0.25;
    end
    return
end

% Checking if base line is stable
d = sr;     % in 1s segments
dd = 1;
mn = [];
time = (0:lenu-1) / sr;
while dd < time(end)
    mn(end+1) = mean(unit(dd:dd+d));
    dd = dd + d;
end
mn(end+1) = mean(unit(dd:end));
if std(mn) > 0.1
    str = 'Instable base line.';
    disp(str)
end

% Segmenting the unit
MS = [];
MI = [];
T = [];
ind2 = 0;
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
    while (npc > nsup || npc < ninf) && next <= 200
        if npc > ninf
            m = m - (m - z) / 2;
        else
            m = m + (m - z) / 2;
        end
        f = find(unitsegment<m);
        npc = length(f) / length(unitsegment);
        next = next + 1;
        if next == 201
            str = 'Maximum number of iterations reached in noise limit search.';
            disp(str)
        end
    end
    MI(end+1) = m;
    
    % Threshold
    T(end+1) = MI(end) + (MS(end) - MI(end)) / 3;
end