function [gr icp err] = linefit2(x,y)
%LINEFIT2   Fits line on points using exact distances for minimalization.
%   [GR ICP ERR] = LINEFIT2(X,Y) fits line on the point-pairs in X and Y and
%   returns the gradient (GR) and the intercept (ICP) of the fitted line as
%   well as the L2-error (ERR) of the fit. Fitting is calculated by 
%   minimizing the exact euclidean distances of the points from the line
%   (in cotrast to POLYFIT, which measures distance along the y-axis).
%
%   LINEFIT2 chooses the shorter between x-wise and y-wise distance for
%   perpendicular distance calculation to avoid stability issues for
%   near-horizontal lines.
%
%   Reference:
%   Sardelis, D. and Valahas, T. "Least Squares Fitting-Perpendicular
%   Offsets
%   (http://mathworld.wolfram.com/LeastSquaresFittingPerpendicularOffsets.html)
%
%   See also POLYFIT and CBSLOPE.

% Fit line
n = length(y);
sy = mean(y);
sx = mean(x);
syy = sum(y.^2);
sxx = sum(x.^2);
B = ((syy - n * sy^2) - (sxx - n * sx^2)) / (2 * ((n * sx * sy) - sum(x.*y)));
b1 = - B + sqrt(B^2+1);
b2 = - B - sqrt(B^2+1);
a1 = sy - b1 * sx;
a2 = sy - b2 * sx;

% Error - deal with stability issues around 0 gradient by choosing the
% better from x-wise and y-wise distance
xwisedist1 = abs(x-(y-a1)/b1);
ywisedist1 = abs(y-(a1+x.*b1));
if ywisedist1(1) < xwisedist1(1)
    Rsquare1 = sum(ywisedist1.^2)/(1+b1.^2);
else
    Rsquare1 = sum((xwisedist1*b1).^2)/(1+b1^2);
end
xwisedist2 = abs(x-(y-a2)/b2);
ywisedist2 = abs(y-(a2+x.*b2));
if ywisedist2(1) < xwisedist2(1)
    Rsquare2 = sum(ywisedist2.^2)/(1+b2.^2);
else
    Rsquare2 = sum((xwisedist2*b2).^2)/(1+b2^2);
end
if Rsquare1 < Rsquare2
    b = b1;
    a = a1;
    Rsquare = Rsquare1;
else
    b = b2;
    a = a2;
    Rsquare = Rsquare2;
end
R = sqrt(Rsquare);

% Output parameters
gr = b;
icp = a;
err = R;