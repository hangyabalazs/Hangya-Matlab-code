function R = lincor(X,Y)
%LINCOR   Linear Correlation.
%   R = LINCOR(X,Y) calculates correlation coefficient between X and Y.
%
%   See also CORRCOEF.

pR = corrcoef(X,Y);
R = pR(2);