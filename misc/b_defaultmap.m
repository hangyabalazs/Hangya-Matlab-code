function M = b_defaultmap(n)
%DEFAULTMAP Creates variable length 'default' colormap.
%   DEFAULTMAP(N) creates 2^N-by-three 'default' colormap. N should be at least 3.
%
%   See also COLORMAP.

%Input arguments check
error(nargchk(1,1,nargin));

m = 2^n;
M = zeros(m,3);

%1st section
c = m/8;
if fix(m) ~= m
    error('N should be at least 3.')
end
for q = 1:c
    M(q,3) = 0.5 + q * (0.5 / c);
end

%2nd section
for w = 1:(2*c)
    M(w+c,2) = w * (0.5 / c);
    M(w+c,3) = 1;
end

%3rd section
for e = 1:(2*c)
    M(e+3*c,1) = e * (0.5 / c);
    M(e+3*c,2) = 1;
    M(e+3*c,3) = 1 - (e - 1) * (0.5 / c);
end

%4th section
for r = 1:(2*c)
    M(r+5*c,1) = 1;
    M(r+5*c,2) = 1 - (r - 1) * (0.5 / c);
end

%5th section
for t = 1:c
    M(t+7*c,1) = 1 - (t - 1) * (0.5 / c);
end