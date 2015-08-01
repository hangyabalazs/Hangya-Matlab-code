function [gr icp err] = linefit(x,y)
%LINEFIT   Fits line on points using exact distances for minimalization.
%   [GR ICP ERR] = LINEFIT(X,Y) fits line on the point-pairs in X and Y and
%   returns the gradient (GR) and the intercept (ICP) of the fitted line as
%   well as the L2-error (ERR) of the fit. Fitting is calculated by 
%   minimizing the exact euclidean distances of the points from the line
%   (in cotrast to POLYFIT, which measures distance along the y-axis).
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
Rsquare1 = sum((y-(a1+x.*b1)).^2)/(1+b1.^2);
Rsquare2 = sum((y-(a2+x.*b2)).^2)/(1+b2.^2);
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