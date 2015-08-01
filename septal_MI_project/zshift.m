function [Z,p] = zshift(eeg,vdisc)
%ZSHIFT    Rayleigh's Z-statistic for shifted unit.
%   [Z,P] = ZSHIFT(EEG,UNIT) returns an array of Z-values calculated with
%   Rayleigh's test for uniformity of phase values of shifted units. The
%   offset is varied from -1 to +1 second. P contains the corresponding
%   significance values.

%   Balazs Hangya
%   Institute of Experimental Medicine
%   Szigony street 43, Budapest 1083, Hungary
%   hangyab@koki.hu

% Hilbert transformation of the eeg
ahee = angle(hilbert(eeg));

% Shift unit
T = -10000:10:10000;
[Z p] = deal(nan(size(T)));
cntr = 0;
for t = T
    cntr = cntr + 1;
    vt = vdisc(vdisc+t>0&vdisc+t<length(eeg))+t;  % shift unit
    bang = ahee(vt);

% Mean resultant length
    n = length(bang);
    ftm = sum(exp(1).^(1i*bang)) / n;    % first trigonometric moment
    mrl = abs(ftm);     % mean resultant length
    z = n * (mrl ^ 2);  % Rayleigh's Z statistic
    Z(cntr) = z;
    p(cntr) = exp(1) ^ (-1 * z) * (1 + (2 * z - z ^ 2) / ...
        (4 * n) - (24 * z - 132 * z ^ 2 + 76 * z ^ 3 - 9 * z ^ 4) / (288 * n ^ 2));
end

% Plot result
figure
plot(T,Z)

sl = 0.05;  % level of significance
siglev = -1 * log(sl);
hold on
line(xlim,[siglev siglev],'Color','k','LineStyle','--');

zmaxloc = find(Z==max(Z));
x_lim = xlim;
y_lim = ylim;
text(x_lim(1)+2*(x_lim(2)-x_lim(1))/3,y_lim(1)+2*(y_lim(2)-y_lim(1))/3,num2str(zmaxloc(1)));