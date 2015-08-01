function df = circdiff(P2,P1,degrad)
%CIRCDIFF   Circular difference.
%   DF = CIRCDIFF(P',P2,S) calculates P2 - P1 on the circle. The value of
%   the string input argument S can be 'deg' or 'rad', determining that the
%   numerical inputs and the output is measured in degrees or in radians.
%
%   See also CIRCULAR_MEAN.

% Conversion between degrees and radians
if isequal(degrad,'rad')
    P1 = P1 * 180 / pi;
    P2 = P2 * 180 / pi;
end

% Difference
df1 = P2 - P1;
if df1 < 0
    df2 = 360 + df1;
else
    df2 = df1 - 360;
end
mdf = min(abs(df1),abs(df2));
if abs(df1) == mdf
    df = df1;
elseif abs(df2) == mdf
    df = df2;
elseif isnan(mdf)
    df = NaN;
else
    error('Technical error 197')
end

% Conversion between degrees and radians
if isequal(degrad,'rad')
    df = df / 180 * pi;
end