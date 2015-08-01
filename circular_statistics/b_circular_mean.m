function av = b_circular_mean(angs,dim)
%CIRCULAR_MEAN   Circular mean.
%   AV = CIRCULAR_MEAN(ANGS) calculates mean of angle values given in radians (ANGS).
%
%   CIRCULAR_MEAN(ANGS,DIM) takes the mean along the dimension DIM of ANGS.
%
%   See also ANGLE, MEAN, MVL, CIRCULAR_STD and CIRCULAR_VAR.

% Input argument check
error(nargchk(1,2,nargin));
mx = max(angs);
mn = min(angs);
x = max(abs(mx),abs(mn));
if x > pi
    error('Input argument must contain radian values.')
end

% Transpose if DIM is 2
trp = 0;
if nargin > 1
    switch dim
    case 1
    case 2
        angs = angs';
        trp = 1;
    otherwise
        error('Unsupported value for dimension.')
    end
end

% Circular mean calculation
S = sin(angs);
C = cos(angs);
sumS = sum(S);
sumC = sum(C);
av = atan(sumS./sumC);

% Transpose back if DIM is 2
if trp
    av = av';
end