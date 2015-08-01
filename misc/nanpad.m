function Xnp = nanpad(X,rd)
%NANPAD   Pad matrix with NaNs.
%   Y = NANPAD(X,N) pads X with N rows and columns of NaNs on both sides,
%   top and bottom.

% Padding
[s1 s2] = size(X);
nans1 = zeros(2*rd+s1,rd) / 0;
nans2 = zeros(rd,s2) / 0;
Xnp = [nans1 [nans2; X; nans2] nans1];