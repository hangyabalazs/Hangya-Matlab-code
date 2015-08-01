function [Rxy_c_z pval] = lincirc_parcor(X,Y,Z)
%LINCIRC_PARCOR   Linear-circular Partial Correlation.
%   R = LINCIRC_PARCOR(X,Y,Z) calculates R(X,Y|Z) partial correlation where
%   X and Z are linear datasets and Y is circular sample.
%
%   See also LINCIRC_CORR, PARCOR and PARTIALCORR.

% Linear and circular correlations
Rxy = lincirc_corr(X,Y);
pR = corrcoef(X,Z);
Rxz = pR(2);
Ryz = lincirc_corr(Z,Y);

% Partial correlation
Rxy_c_z = (Rxy - Rxz * Ryz) / (sqrt(1-Rxz^2) * (sqrt(1-Ryz^2)));

% Test for non-zero correlation
[n dz] = size(Z);
df = max(n - dz - 2,0); % degrees of freedom
t = sign(Rxy_c_z) .* Inf;
k = (abs(Rxy_c_z) < 1);
t(k) = Rxy_c_z(k) ./ sqrt(1-Rxy_c_z(k).^2);
t = sqrt(df).*t;
tail = 'b';
switch tail
    case 'b' % 'both or 'ne'
        pval = 2*tcdf(-abs(t),df);
    case 'r' % 'right' or 'gt'
        pval = tcdf(-t,df);
    case 'l' % 'left or 'lt'
        pval = tcdf(t,df);
end