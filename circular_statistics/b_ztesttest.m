function z = b_ztesttest

% Angle
x = [0:0.01:16*pi];
ang = angle(cos(x)+i*sin(x));

% Random process
% r = rand(1,1000);
% cs = cumsum(r);
% ind = find(cs<max(x));
% ind = ind(end);
% cs = cs(1:ind);
cs = linspace(1,max(x),1000);

% Rayleigh's Z-test
bang = ang(ceil(cs));

% Mean resultant length
n = length(bang);
ftm = sum(exp(1).^(i*bang)) / n;    % first trigonometric moment
mrl = abs(ftm);     % mean resultant length
z = n * (mrl ^ 2);  % Rayleigh's Z statistic
p = exp(1) ^ (-1 * z) * (1 + (2 * z - z ^ 2) / ...
    (4 * n) - (24 * z - 132 * z ^ 2 + 76 * z ^ 3 - 9 * z ^ 4) / (288 * n ^ 2));