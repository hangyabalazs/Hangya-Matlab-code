function R = mvltest
%MVLTEST   Testing continuous extension of mean resultant length.


phi = 0:pi/16:2*pi;
phi = repmat(phi,1,8);
r = 2 + (randn(size(phi)) * 0.15);

R = 1 / sum(r) * abs(sum(r.*exp(i*phi)));