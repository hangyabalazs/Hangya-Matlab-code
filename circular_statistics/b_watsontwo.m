function [u2 p] = b_watsontwo(x,y)
%WATSONTWO    Watson's test for homogeneity on two samples of circular data.
%   [U2,P] = WATSONTWO(X,Y) performs test for homogeneity on vectors of 
%   circular data measured in radians (X and Y). Test statistic (U2) and
%   p-value (P) is returned.
%
%   Critical values for the test statistic are obtained using the
%   asymptotic distribution of the test statistic.  It is recommended to
%   use the obtained critical values and ranges for p-values only
%   for combined sample sizes in excess of 17.
%
%   See also WATSON and RAO.

% Input argument check
error(nargchk(2,2,nargin))

% Sample size
n1 = length(x);
n2 = length(y);
n = n1 + n2;
if n < 18
    error('Sample size is too small.')
end

% Test
x = [sort(mod(x,2*pi)), repmat(1,n1,1)];
y = [sort(mod(y,2*pi)), repmat(2,n2,1)];
xx = [x; y];
[S,rank] = sort(xx(:,1));
xx = [xx(rank,:) (1:n)'];
a = [1:n];
b = [1:n];
for i = 1:n
    a(i) = sum(xx(1:i,2)==1);
    b(i) = sum(xx(1:i,2)==2);
end
d = b / n2 - a / n1;
dbar = mean(d);
u2 = (n1 * n2) / n^2 * sum((d-dbar).^2);
crits = [99, 0.385, 0.268, 0.187, 0.152];
if u2 > 0.385
    p(1) = 0;
    p(2) = 0.001;
elseif u2 > 0.268
    p(1) = 0.001;
    p(2) = 0.01;
elseif u2 > 0.187
    p(1) = 0.01;
    p(2) = 0.05;
elseif u2 > 0.152
	p(1) = 0.05;
    p(2) = 0.10;
else
    p(1) = 0.10;
    p(2) = 1;
end