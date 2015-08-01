function Rxy_c_z = parcor(X,Y,Z)
%PARCOR   Partial Correlation.
%   R = PARCOR(X,Y,Z) calculates the partial correlation R(X,Y|Z) similar
%   to PARTIALCORR (Statistics Toolbox function).
%
%   See also PARTIALCORR

% Linear correlation
pR = corrcoef(X,Y);
Rxy = pR(2);
pR = corrcoef(X,Z);
Rxz = pR(2);
pR = corrcoef(Y,Z);
Ryz = pR(2);

% Partial correlation
Rxy_c_z = (Rxy - Rxz * Ryz) / (sqrt(1-Rxz^2) * (sqrt(1-Ryz^2)));