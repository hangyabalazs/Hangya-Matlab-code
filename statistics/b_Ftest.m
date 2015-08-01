function [stat,h,p] = b_Ftest(d1,d2,alpha)
%FTEST   F-test.
%   [STAT,H,P] = FTEST(D1,D2) returns the F-stat (STAT) and the
%   corresponding p (P) value for D1 and D2. H equals to 0 if
%   nullhypothesis of equal standard deviation holds at 5% significance
%   level, H equals to 1 otherwise.
%
%   FTEST(D1,D2,ALPHA) uses ALPHA signifance level for hypothesis testing.
%
%   See also TTEST and TTEST2.

% Input argumnet check
error(nargchk(2,3,nargin))
if nargin == 2
    alpha = 0.05;
end

% F-test
s1 = max(var(d1),var(d2));
s2 = min(var(d1),var(d2));
stat = s1 / s2;
v1 = length(d1) - 1;
v2 = length(d2) - 1;
if isequal(s1,var(d1))
    p = 2 * (1 - fcdf(stat,v1,v2));
elseif isequal(s1,var(d2))
    p = 2 * (1 - fcdf(stat,v2,v1));
end
if p > alpha
    h = 0;
else
    h = 1;
end