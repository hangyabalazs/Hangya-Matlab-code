function [Z,p] = b_zshift2
%ZSHIFT2    Rayleigh's Z-statistic for shifted unit.
%   [Z,P] = ZSHIFT2 returns an array of Z-values calculated with Rayleigh's
%   test for uniformity of phase values of shifted units. The offset is
%   varied from -1 to +1 second. P contains the correspon

% Import and discrimination
b_in3
b_disc
b_var2ws('all','caller')

% Hilbert transformation of the eeg
ahee = angle(hilbert(eeg));
aheedeg = ahee * (180 / pi);

% Shift unit
Z = [];
p = [];
T = [-10000:10:10000];
for t = T
    vt = vdisc(find(vdisc+t>0&vdisc+t<length(eeg)))+t;  % shift unit
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
figure
plot(T,Z)

sl = 0.05;  % level of significance
% siglev = exp(1) ^ (-1 * sl) * (1 + (2 * sl - sl ^ 2) / ...
%     (4 * n) - (24 * sl - 132 * sl ^ 2 + 76 * sl ^ 3 - 9 * sl ^ 4) / (288 * n ^ 2));
siglev = -1 * log(sl);
hold on
line(xlim,[siglev siglev],'Color','k','LineStyle','--');

zmaxloc = find(Z==max(Z));
x_lim = xlim;
y_lim = ylim;
text(x_lim(1)+2*(x_lim(2)-x_lim(1))/3,y_lim(1)+2*(y_lim(2)-y_lim(1))/3,num2str(zmaxloc(1)));