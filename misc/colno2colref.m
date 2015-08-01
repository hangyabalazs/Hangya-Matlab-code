function colref = colno2colref(colno)
%COLNO2COLREF   Converts column index to Excel column reference.
%   COLREF = COLNO2COLREF(COLNO) converts column index (COLNO) to
%   corresponding Excel column reference (COLREF).
%
%   See also XLSREAD and XLSWRITE.

% Convert column index to Excel column reference
S = ['A' 'B' 'C' 'D' 'E' 'F' 'G' 'H' 'I' 'J' 'K' 'L' 'M' 'N' 'O' 'P' 'Q' 'R' 'S' 'T' 'U' 'V' 'X' 'Y' 'Z'];
n = length(S);
r = mod(colno,n);
if r == 0
    r = n;
end
q = floor((colno-1)/n);
if q == 0
    colref = S(r);
else
    colref = [S(q) S(r)];
end