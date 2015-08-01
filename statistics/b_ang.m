function Phi = b_ang(Z1,Z2)
%ANG    Computes the angle of two vectors.
%   Phi = ANG(Z1,Z2) computes the angle of vetors Z1 and Z2 in the usual
%   euclidian space.

% Input arguments check
error(nargchk(2,2,nargin));

% Computing the angle
N1 = norm(Z1,2);
N2 = norm(Z2,2);
ScalarProduct = b_scalp(Z1,Z2);
CosPhi = ScalarProduct / (N1 * N2);
Phi = acos(CosPhi);