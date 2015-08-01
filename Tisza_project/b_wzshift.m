function [zmaxloc_t,ang_t,mvl_t zmaxloc_c,ang_c,mvl_c] = b_wzshift(tI,cI,vdtI,vdcI)
%WZSHIFT   Z-shift method for Tisza data.
%   WZSHIFT implements z-shift method for Tisza runoff data. See ZSHIFT and
%   ZSHIFTRUN3 for detailed description on the z-shift method.
%   Input arguments:
%       tI: Tisza runoff data
%       cI: Szamos runoff data
%       vdtI: discriminated tI
%       vdcI: discriminated cI
%
%   Output arguments:
%       zmaxloc_t: z-shift value for Tisza events
%       ang_t: mean angle for Tisza phases
%       mvl_t: mean vector length for Tisza phases
%       zmaxloc_c: z-shift value for Szamos events
%       ang_c: mean angle for Szamos phases
%       mvl_c: mean vector length for Szamos phases
%
%   See also ZSHIFT, ZSHIFTRUN2 and WZSHIFT2.

[zmaxloc_t,ang_t,mvl_t] = zshift(vdtI,cI);
[zmaxloc_c,ang_c,mvl_c] = zshift(vdcI,tI);

% -------------------------------------------------------------------------
function [zmaxloc,ang,mvl] = zshift(vdisc,eeg)

% Filtering
flt = fir1(512,[26/(365/2) 52/(365/2)]);        % 1-2 weeks period time
feeg = filtfilt(flt,1,eeg);

% Hilbert transformation of the eeg
ahee = angle(hilbert(feeg));
aheedeg = ahee * (180 / pi);

% Phase angles
bahee = ahee(vdisc);
n = length(bahee);
ftm = sum(exp(1).^(i*bahee)) / n;    % first trigonometric moment
ang = angle(ftm);   % mean angle
mvl = abs(ftm);     % mean resultant length

% Shift unit
Z = [];
p = [];
T = [-30:1:30];
for t = T
    vt = vdisc(find(vdisc+t>0&vdisc+t<length(eeg)))+t;  % shift unit
    if isempty(vt)      % Skip segment, if all spikes fall out of range after unit shift
        zmaxloc = NaN;
        figure(H);
        plot(0,0)
        text(0,0,'All spikes out of range after unit shift.','Color','red',...
            'HorizontalAlignment','center');
        return
    end
    bang = ahee(vt);

% Mean resultant length
    n = length(bang);
    ftm = sum(exp(1).^(i*bang)) / n;    % first trigonometric moment
    mrl = abs(ftm);     % mean resultant length
    z = n * (mrl ^ 2);  % Rayleigh's Z statistic
    Z(end+1) = z;
    p(end+1) = exp(1) ^ (-1 * z) * (1 + (2 * z - z ^ 2) / ...
        (4 * n) - (24 * z - 132 * z ^ 2 + 76 * z ^ 3 - 9 * z ^ 4) / (288 * n ^ 2));
end

% Plot result
figure;
plot(T,Z)

sl = 0.05;  % level of significance
siglev = -1 * log(sl);
hold on
line(xlim,[siglev siglev],'Color','k','LineStyle','--');

zmxlc = find(Z==max(Z));
x_lim = xlim;
y_lim = ylim;
text(x_lim(1)+2*(x_lim(2)-x_lim(1))/3,y_lim(1)+2*(y_lim(2)-y_lim(1))/3,num2str(T(zmxlc(1))));
hold off

zmaxloc = T(zmxlc);