function [tI cI] = b_wread3
%WREAD3   Import Tisza data.
%
%   See also WDISC and IN4.

% Read
fname = 'd:\_analysis\matlab_data\Tisza\TBECSQ_2.txt';
[year month day hour minute y] = textread(fname,'%f %f %f %f %f %f');
data = [year month day hour minute y];

month2 = [0 31 59 90 120 151 181 212 243 273 304 334];

x = (year - 2001) *  365 * 24 * 60 + month2(month)' * 24 * 60 + (day - 1) * 24 * 60 ...
    + hour * 60 + minute;
[x(1) x(end)]
data(1,:)
xi = [0:15:x(end)];
yi = interp1(x,y,xi);

save('d:\_analysis\matlab_data\Tisza\tbecsq_interp.mat','xi','yi')


% x2 = x / 60;
% x3 = x2 / 24;
% figure
% plot(x)

% figure
% b_lombper(x(1:2000),y(1:2000));
% figure
% b_lombper(x2(1:2000),y(1:2000));
% figure
% b_lombper(x3(1:2000),y(1:2000));