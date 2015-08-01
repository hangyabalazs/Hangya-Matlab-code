function b_retimefig(time);
%RESCALEFIG   Changes ticklabeles from 'matrix index' to 'time'.
%   RETIMEFIG(TIME) requires an input parameter TIME containing the correspondig
%   time value for each matrix index.
%
%   See also RESCALEFIG and WAVELET.

a0 = get(gca,'XTickLabel');
a00 = str2num(a0);
iszero = 0;
if a00(1) == 0
    a00(1) = 1;
    iszero = 1;
end
a1 = time(a00);
if iszero
    a1(1) = 0;
end
set(gca,'XTickLabel',num2str(a1'))