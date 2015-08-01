function r = mod2(a,m)
%MOD2   Modulus.
%   R = MOD2(A,M) operates as MOD but returns M instead of 0.
%
%   See also MOD.

r = mod(a,m);
r(r==0) = m;