function corrmap = czcorrmap(rm1,rm2)
%CZCORRMAP   Correlation map.
%   C = CZCORRMAP(R1,R2) returns the correlation map of R1 and R2 in C.
%   Linear correlation is calculated between 10x10 parts of the input
%   matrices.
%
%   See also CZPLACE.

% Input argument check
error(nargchk(2,2,nargin))

% Calculate correlation map
s1 = size(rm1,1);
s2 = size(rm1,2);
corrmap = zeros(size(rm1));
for x = 1:s1
    for y = 1:s2
        if ~isnan(rm1(x,y)) && ~isnan(rm2(x,y))
            xind1 = max(x-5,1);
            xind2 = min(x+5,s1);
            yind1 = max(y-5,1);
            yind2 = min(y+5,s2);
            c1 = rm1(xind1:xind2,yind1:yind2);
            c2 = rm2(xind1:xind2,yind1:yind2);
            c11 = c1(~isnan(c1));
            c21 = c2(~isnan(c2));
            C = corrcoef(c11,c21);
            corrmap(x,y) = C(2);
        else
            corrmap(x,y) = NaN;
        end
    end
end

% Plot
figure
pcolor(corrmap)