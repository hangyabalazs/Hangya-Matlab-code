function v = b_circular_var(angs,dim)
%CIRCULAR_VAR   Circular variance.
%   V = CIRCULAR_STD(ANGS) calculates variance for angle values given in radians (ANGS).
%
%   CIRCULAR_STD(ANGS,DIM) takes the variance along the dimension DIM of ANGS.
%
%   See also ANGLE, STD, CIRCULAR_MEAN, CIRCULAR_STD and MVL.

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

% Circular STD calculation
S = sin(angs);
C = cos(angs);
sumS = sum(S);
sumC = sum(C);
n = length(angs);
arg = sumS .^2 + sumC .^ 2;
len = sqrt(arg) / n;
v = 1 - len;

% Transpose back if DIM is 2
if trp
    v = v';
end