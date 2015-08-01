function R = lincirc_corr(x,y)
%LINCIRC_CORR   Linear-circular correlation.
%   LINCIRC_CORR(X,Y) calculates linear-circular correlation of X (linear)
%   and Y (circular) variables.
%
%   Reference: Mardia KV, Jupp PE (2000) Directional statistics. Wiley, 
%   Chichester, UK, pp. 245-246.

rxc = lcorr(x,cos(y));
rxs = lcorr(x,sin(y));
rcs = lcorr(cos(y),sin(y));
Rsquare = (rxc^2 + rxs^2 - 2 * rxc * rxs * rcs) / (1 - rcs^2);
R = sqrt(Rsquare);

% -------------------------------------------------------------------------
function R = lcorr(x,y)

pR = corrcoef(y',x');
R = pR(2);