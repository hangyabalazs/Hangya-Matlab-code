function dev = b_circular_std(angs,dim)
%CIRCULAR_STD   Circular standard deviation.
%   DEV = CIRCULAR_STD(ANGS) calculates standard deviation for angle values given
%   in radians (ANGS).
%
%   CIRCULAR_STD(ANGS,DIM) takes the standard deviation along the dimension DIM of 
%   ANGS.
%
%   See also ANGLE, STD, CIRCULAR_MEAN, CIRCULAR_VAR and MVL.

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
arg1 = sumS .^2 + sumC .^ 2;
arg2 = sqrt(arg1) / n;
arg3 = -2 * log(arg2);
dev = sqrt(arg3);


% Transpose back if DIM is 2
if trp
    dev = dev';
end