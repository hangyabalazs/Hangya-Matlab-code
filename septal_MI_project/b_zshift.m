function [Z,p,Z_minus,p_minus,Z_plus,p_plus] = b_zshift

% Import and discrimination
b_in3
b_disc
b_var2ws('all','caller')

% Shift unit
vdisc_minus = vdisc(find(vdisc-10000>0))-10000;
vdisc_plus = vdisc(find(vdisc+10000<length(eeg)))+10000;

% Hilbert transformation of the eeg
ahee = angle(hilbert(eeg));
aheedeg = ahee * (180 / pi);
bang = ahee(vdisc);
bang_minus = ahee(vdisc_minus);
bang_plus = ahee(vdisc_plus);

% Mean resultant length
n = length(bang);
ftm = sum(exp(1).^(i*bang)) / n;    % first trigonometric moment
mrl = abs(ftm);     % mean resultant length
Z = n * (mrl ^ 2);  % Rayleigh's Z statistic
p = exp(1) ^ (-1 * Z) * (1 + (2 * Z - Z ^ 2) / ...
    (4 * n) - (24 * Z - 132 * Z ^ 2 + 76 * Z ^ 3 - 9 * Z ^ 4) / (288 * n ^ 2));

n = length(bang_minus);
ftm = sum(exp(1).^(i*bang_minus)) / n;    % first trigonometric moment
mrl = abs(ftm);     % mean resultant length
Z_minus = n * (mrl ^ 2);  % Rayleigh's Z statistic
p_minus = exp(1) ^ (-1 * Z_minus) * (1 + (2 * Z - Z ^ 2) / ...
    (4 * n) - (24 * Z - 132 * Z ^ 2 + 76 * Z ^ 3 - 9 * Z ^ 4) / (288 * n ^ 2));

n = length(bang_plus);
ftm = sum(exp(1).^(i*bang_plus)) / n;    % first trigonometric moment
mrl = abs(ftm);     % mean resultant length
Z_plus = n * (mrl ^ 2);  % Rayleigh's Z statistic
p_plus = exp(1) ^ (-1 * Z_plus) * (1 + (2 * Z - Z ^ 2) / ...
    (4 * n) - (24 * Z - 132 * Z ^ 2 + 76 * Z ^ 3 - 9 * Z ^ 4) / (288 * n ^ 2));