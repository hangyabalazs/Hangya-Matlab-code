function C = b_ifelse(x,yes,no)
%IFELSE    Conditional mixing.
%   IFELSE(X,YES,NO) mixes the elements of YES and NO in the following way:
%   it takes the element of YES if corresponding element of logical input X
%   is true, otherwise it takes the corresponding element of NO. It returns
%   an X-size matrix.
%
%   Note that YES and NO should be either X-size, or 1x1!
%
%   See also A1INV.

% Input argument check
error(nargchk(3,3,nargin)) 
if ~isequal(size(x),size(yes))
    if ~isequal(size(yes),[1 1])
        error('Input argumnets X and YES must be of equal size or YES should be 1-by-1.')
    else
        yes = repmat(yes,size(x));
    end
end
if ~isequal(size(x),size(no))
    if ~isequal(size(no),[1 1])
        error('Input argumnets X and NO must be of equal size or NO should be 1-by-1.')
    else
        no = repmat(no,size(x));
    end
end
    
% Mixing
A = x .* yes;
B = ~x .* no;
C = A + B;