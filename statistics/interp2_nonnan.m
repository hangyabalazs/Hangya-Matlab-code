function Y2 = interp2_nonnan(X,nr)
%INTERP2_NONNAN   2D interpolation.
%   INTERP2_NONNAN(X,N) calls INTERP2 with input arguments X and N (see 
%   INTERP2 for details). It handles NaNs in X: interpolates at NaN 
%   locations and restores NaNs at the end.
%
%   See also INTERP2.

% Interpolate NaNs
X2 = X;
siz = size(X);
while any(isnan(X2(:)))
    fn = find(isnan(X2));
    for k = 1:length(fn)
        [x y] = ind2sub(siz,fn(k));
        inx1 = max(1,x-1);
        iny1 = max(1,y-1);
        inx2 = min(x+1,siz(1));
        iny2 = min(y+1,siz(2));
        subX = X2(inx1:inx2,iny1:iny2);
        if ~all(all(isnan(subX)))
            X2(fn(k)) = b_mean_nonnan(subX);
        end
    end
end

% 2D interpolation
Y1 = interp2(X2,nr);

% Restore NaNs
Y2 = Y1;
fn = find(isnan(X));
for k = 1:length(fn)
    [x y] = ind2sub(siz,fn(k));
    newx2 = x * 2 ^ nr - 2 ^ nr + 1;
    newy2 = y * 2 ^ nr - 2 ^ nr + 1;
    newx1 = min(newx2+2^nr-1,size(Y2,1));
    newy1 = min(newy2+2^nr-1,size(Y2,2));
    Y2(newx2:newx1,newy2:newy1) = NaN;
end