%Test program for linear-circular partial correlation.

lim = 1;
X = (-lim:0.001:lim)';        % phase;
X = (2*rand(1000,1)-1)*3.1
X = rand(1000,1)*2*3.1
Y = 2*X + rand(length(X),1); % RT
Z = 2*X + rand(length(X),1); % amp

disp(' ')
disp(' ')
disp(' ')

% Linear correlation
pR = corrcoef(X,Y);
Rxy = pR(2);
Rxy

% Circular correlation
% cRxy = lincirc_corr(Y',X')
cRxy = lincirc_corr2(Y,X)
% pR = corrcoef([sin(2*atan(Y)) cos(2*atan(Y))],[sin(X) cos(X)]);
% cRxy = pR(2)
% pR = corrcoef([Y Y],[sin(X) cos(X)]);
% cRxy = pR(2)

% Partial correlation (linear)
Rxy_c_z = parcor(X,Y,Z)
% Rxy_c_z / Rxy

% Partial correlation (circular)
% cRxy_c_z = lincirc_parcorr(Y',X',Z')
cRxy_c_z = lincirc_parcorr5(Y',X',Z')