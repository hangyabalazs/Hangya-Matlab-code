function Ydp = smooth2_nonnan(X,gk)
%SMOOTH2_NONNAN   2D smoothing.
%   Y = SMOOTH2_NONNAN(X,K) smoothes X using K as convolution Kernel and
%   returnes the smoothed matrix in Y. It pads X with NaNs and ignores NaNs
%   when convolving (weights are renormalized).
%
%   See also SMOOTH2 and GAUSSKERNEL.

% Input argument check
error(nargchk(2,2,nargin))

% Padding
[s1 s2] = size(gk);
rd = floor(max(s1,s2)/2);    % radius
Xnp = nanpad(X,rd);
[snp1 snp2] = size(Xnp);

% Smoothing
Y = zeros(size(Xnp));
for x = 1:snp1
    for y = 1:snp2
        if ~isnan(Xnp(x,y))
            xind1 = x - rd;
            xind2 = x + rd;
            yind1 = y - rd;
            yind2 = y + rd;
            C = Xnp(xind1:xind2,yind1:yind2);
            ngk = gk .* double(~isnan(C));
            ngk = ngk / sum(sum(ngk));
            prd = C .* ngk;
            Y(x,y) = b_sum_nonnan(prd(:));
        else
            Y(x,y) = NaN;
        end
    end
end

% Remove the effect of padding
Ydp = depad(Y,rd);



% -------------------------------------------------------------------------
function Xnp = nanpad(X,rd)

[s1 s2] = size(X);
nans1 = zeros(2*rd+s1,rd) / 0;
nans2 = zeros(rd,s2) / 0;
Xnp = [nans1 [nans2; X; nans2] nans1];



% -------------------------------------------------------------------------
function Ydp = depad(Y,rd)

Ydp = Y(rd+1:end-rd,rd+1:end-rd);