function [speed_vec , time_vec] = get_speed(filepath)
% FUNCTION [SPEED_VEC, TIME_VEC] = GET_SPEED(FILEPATH)

if nargin<1,
    help get_speed
    loadcb
    allcells = listtag('cells');
    filepath = [cellid2fnames(allcells(100),'Sess')  filesep 'VT1.mat'];
end

% PARAMETER INITIALIZATION
bns = 8;                            % number of pixels to bin/jump to allow.
numbins = 20;                       % number of bins for speed trajectory.
linearX = 0;                        % Consider only changes in X position.
linearY = 0;                        % Consider only changes in Y position.
ifsmooth = 1;                       % Do you want to smooth the speed.
span = 5;                           % window span for moving window for speed.
win_margin = [0 0];             % margin around the first and second event.
xpix = 1;                               % pixels / cm
ypix = 1;                               % pixels/cm

try
    load(filepath);
catch
    error('Could not locate position file for cellid')
end
x_pos=ExtractedX;
y_pos=ExtractedY;

if linearX ==1,
    y_pos = ones(size(ExtractedY));
elseif linearY == 1,
    x_pos = ones(size(ExtractedX));
else
end
x_pos_ts=TimeStamps*1e-6;

% Position vector correction
dx = diff(x_pos);       % leave out the outliers
fdx = find(abs(dx)>bns&abs([dx(2:end) 0])>bns&abs(dx-[dx(2:end) 0])>2*bns);
x_pos(fdx+1) = [];
y_pos(fdx+1) = [];
x_pos_ts(fdx+1) = [];
dy = diff(y_pos);
fdy = find(abs(dy)>bns&abs([dy(2:end) 0])>bns&abs(dy-[dy(2:end) 0])>2*bns);
x_pos(fdy+1) = [];
y_pos(fdy+1) = [];
x_pos_ts(fdy+1) = [];

inxs = x_pos > 0 & y_pos > 0;   % leave out (0;0) points
X = x_pos(inxs);
Y = y_pos(inxs);
TP = x_pos_ts(inxs);

% Interpolate the deleted outlier timepoints.
TP2 = TimeStamps*1e-6;
X = interp1(TP,X,TP2);
Y = interp1(TP,Y,TP2);
TP = TP2;
% Note: dt is not accurate: 0.0204 : 0.0459 (twice sometimes)
dt_vec = diff([0 TP]);
clear ExtractedX ExtractedY TimeStamps TP2

Xdist = diff([0 X])./xpix; % cm
Ydist = diff([0 Y])./ypix; % cm

speed = sqrt(Xdist.^2 + Ydist.^2)./dt_vec;

% calculate speed.
% speed = sqrt(diff([0 X]).^2 + diff([0 Y]).^2);
% speed = sqrt(diff([0 X]).^2 + diff([0 Y]).^2)./dt_vec;
% speed3 = sqrt([0 diff(X)].^2 + [0 diff(Y)].^2 )./dt_vec;

if ifsmooth == 1,
    speed_vec = smooth(speed,span,'moving');
else
    speed_vec = speed;
end
time_vec = TP';

% -------------------------------------------------------------------------
function [X2 S] = smooth(X,str,wn)
%SMOOTH   Smoothing by a moving average.
%   [O S] = SMOOTH(X,M,WN) performs smoothing on X by averaging in a
%   sliding window of WN size. M describes the data type, i.e. 'linear' or 
%   'circular'. Output is returned in O with standard error of the
%   smoothing in S.
%
%   See also CONV.

n = wn;
nn = (n - 1) / 2;
m = length(X);
X2 = zeros(size(X));
S = zeros(size(X));
switch str
    case 'linear'
        for k = 1:m
            ind1 = max(1,k-nn);
            ind2 = min(k+nn,m);
            X2(k) = mean(X(ind1:ind2));
            S(k) = std(X(ind1:ind2));
        end
    case 'circular'
        for k = 1:m
            if k - nn < 1
                X2(k) = mean([X(mod2(k-nn,m):m); X(1:k+nn)]);
                S(k) = std([X(mod2(k-nn,m):m); X(1:k+nn)]) / sqrt(n);
            elseif k + nn > m
                X2(k) = mean([X(k-nn:m); X(1:mod2(k+nn,m))]);
                S(k) = std([X(k-nn:m); X(1:mod2(k+nn,m))]) / sqrt(n);
            else
                X2(k) = mean(X(k-nn:k+nn));
                S(k) = std(X(k-nn:k+nn)) / sqrt(n);
            end
        end
end