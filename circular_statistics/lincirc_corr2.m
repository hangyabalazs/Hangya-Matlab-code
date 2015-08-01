function [R p] = lincirc_corr2(x,y)
%LINCIRC_CORR2   Linear-circular correlation.
%   [R P] = LINCIRC_CORR2(X,Y) calculates linear-circular correlation (R)
%   of X (linear) and Y (circular) variables. Significance level is
%   returned in P.
%
%   Reference: Fisher NI (1993) Statistical analysis of circular data.
%   Cambridge University Press
%
%   See also LINCIRC_CORR.

[b,bint,r,rint,stats] = regress(x,[ones(length(y),1) ...
    sin(y) cos(y)]);      % correlation
R = sqrt(stats(1));
p = stats(3);