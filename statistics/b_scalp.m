function P = b_scalp(Z1,Z2)
%SCALP  Computes scalarproduct of two vectors.
%   P = SCALP(Z1,Z2) computes the usual euclidian scalarproduct
%   of the input vectors Z1 and Z2.

% Input arguments check
error(nargchk(2,2,nargin));
if ~isreal(Z1) | ~isreal(Z2)
    warning('At least one of the input arguments is complex.')
end

% Computing the scalar product
S1 = size(Z1);
S2 = size(Z2);
f1 = find(S1==1);
f2 = find(S2==1);
if isempty(f1) | isempty(f2)
    error('Input arguments must be arrays.')
end
if f1(1) == 2
    Z1 = Z1';
end
if f2(1) == 1
    Z2 = Z2';
end
P = Z1 * Z2;